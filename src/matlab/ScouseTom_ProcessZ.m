function [ Zall ] = ScouseTom_ProcessZ( HDR,TT,ExpSetup,PlotFlag )
%[ Zall ] = ScouseTom_ProcessZ( HDR,TT ) 
%   SCOUSETOM_LOADZ
%   Calculates the contact impedance from a two electrode "contact check"
%   measurement. The impedance is estimated from the voltage on the
%   injection channels. This type of injection is started by the
%   ScouseTom_ContactCheck command for the ScouseTom repository.
%
%   The flow is much the same as ScouseTom_ProcessBV, execpt it it simpler
%   in that the data is small enough to load all at once given the short
%   injection times during contact checks.
%
%   Currently the amplitude and frequency are FIXED, and as such a
%   _log.mat is not stored alongside the EEG datafiles with contact checks.
%
%   As with ProcessBV, the filter settings are found based on the injection
%   frequency and the time between switching electrode pairs. The mean
%   voltage on all channels is then found by demodulating using the Hilbert
%   transform. Given the amplitde of injected current, the contact
%   impedance is estimated from the mean voltage on each electrode when it
%   is used to inject current.
%
%   Bar graph indicating bad channels is produced. Useful during electrode
%   application.
%
%   The output is stored in a matfile Fname-ZCheck.mat
%
%   Inputs:
%   HDR - from ScouseTom_getHDR contains metadata about EEG file
%   TT  - from ScouseTom_TrigProcess, infomation from digital trigger
%       channels
%   ExpSetup -  Structure containing the setup for the ScouseTom system. 
%       Most important variables can be assumed to be default, so only
%       necessary if you ahve altered these in the firmware (optional)
%   PlotFlag[1] - Specify if you wna the bar plot (optional)
%
%   Outputs:
%   Zall - Structure containing the Z structure for every zcheck in the
%   file. 
%
%   Z structure:
%   Z(Chn) - Estimated contact impedance - averaged across both injections on
%       that electrode
%   Zi(Chn,2) - Estimated contact impedance, per injection
%   dZ(Chn) - Difference in voltage, in form of previous UCH system
%   ExpSetup - Input structure for ref
%   TimeNum - timestamp unix format
%   TimeVec - timestamp in matlab timevector format

%% check inputs are ok here

%normally we want to plot
if exist('PlotFlag','var') ==0 || isempty(PlotFlag)
    PlotFlag =1;
end

%set defaults 
if exist('ExpSetup','var') ==0 || isempty(ExpSetup)
    ExpSetup.Elec_num=HDR.NS-1; %number of electrodes
    ExpSetup.Bad_Elec=[]; %bad electrodes assume none
    ExpSetup.Desc='THIS IS A DUMMY EXPSETUP MADE DURING ZCHECK';
    ExpSetup.Amp=141; % fixed in ScouseTom firmware
    ExpSetup.Freq=1000; % fixed in ScouseTom firmware
end

%% See what type of system we are using

%the maximum voltage is different for each system, this *should* be in the
%HDR structure somewhere, but I dont know where it is

switch HDR.TYPE
    case 'BDF' % biosemi file
        MaxV=0.5e6; %500mV range on BioSemi
    case 'BrainVision' %actichamp
        MaxV=0.25e6; %250mV range on ActiChamp
    otherwise
        error('Unknown HDR Type');
end

Fs=HDR.SampleRate;
fname=HDR.FILE.Name;

%% Get variables

Nelec=ExpSetup.Elec_num; %number of electrodes
Nprt=Nelec; %this is the same as N_elec as the protocol lines is equal to number of electrodes for contact checks

%create the contact check protocol - assuming neighbouring pairs
Contact_Protocol=[1:Nelec; circshift(1:Nelec,-1,2)]';
Contact_Protocol=Contact_Protocol(~(any(ismember(Contact_Protocol,ExpSetup.Bad_Elec),2)),:); %remove reference to bad electrodes

%% Amplitude and freq

%these are fixed in the arduino code in Injection.h in ScouseTom repo
Amp=ExpSetup.Amp; % 141 uA
Freq=ExpSetup.Freq; %1 kHz

%scale factor - impedance conversion
ZSF=1/(Amp); %VOLTAGE GIVEN IN uV SO AMP in uA too
Z_Max=MaxV*ZSF; %maximum impedance based on maximum voltage measureable on system

%% Defaults for filtering and contact impedance values

%bandwidth of filter
BW=50;
Z_Rec=1000; %recommended contact impedance

%% run each impedance check separately

% run each
Z_checks_num=length(TT.Contact.InjectionStarts);

if isempty(Z_checks_num)
    warning('No Contact Checks Found');
    return
end

for iZ = 1:Z_checks_num
    
    
    %% read portion of interest
    
    %sread needs integer seconds
    Data_start=floor(TT.Contact.InjectionStarts(iZ)/Fs);
    Data_end=ceil(TT.Contact.InjectionStops(iZ)/Fs);
    Data_length=fix(Data_end-Data_start);
    
    %convert to number of samples
    Data_start_s=Data_start*Fs;
    Data_end_s=Data_end*Fs;
    
    %load the data
    disp(['Reading relevant data for contact check: ', num2str(iZ)]);
    
    %sread MUST take INTEGER seconds as inputs, hence rounding data to
    %nearest second as above
    V=sread(HDR,Data_length,Data_start);
    
    % take relevant switches
    curInjSwitch=TT.Contact.InjectionSwitches{iZ};
    
    %use the second injection as the estimate for the filtering - in case
    %something was wrong with first (warming up or delay starting etc.)
    swidx=2;
    
    %determine the Carrier Frequency, Filter coeffs, and amount of signal
    %to trim from the electrode on the *second* injection
    %find the corresponding filter settings
    [Filt,FilterTrim,Fc]=ScouseTom_FindFilterSettings(HDR,TT.Contact.InjectionSwitches(iZ,:),Contact_Protocol(swidx,1));
    
    %demodulate each segment in turn using hilbert transfrom
    [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,Filt{1}); % filter and demodulate channel
    
    %get just the mean of the voltages during each injection
    [Vmag,PhaseRaw]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,FilterTrim{1},TT.Contact.InjectionSwitches{iZ}-Data_start_s); %process each injection window, adjusting for new start time
    
    %Get Phase
    [Phase]=ScouseTom_data_PhaseEst(PhaseRaw,Contact_Protocol,1);
    
    %% Calculate Z
    CurElecs=unique(Contact_Protocol);
    
    %for each electrode
    for iElec=1:Nelec
        %if it was skipped as a BAD ELEC
        if ismember(iElec,CurElecs)
            %find the voltages when this electrode is injecting
            InjchnV=Vmag(any(Contact_Protocol == iElec,2),iElec);
            %convert into impedances using peak amplitude values and convert uV to V
            Z(iElec,:)=InjchnV*ZSF;
            
            %calc the "legacy" UCH 2 channel one too - the difference
            %between adjacent pairs
            Injchndiff=diff(Vmag(Contact_Protocol(:,1) == iElec,1:2));
            dZ(iElec)=abs(Injchndiff*ZSF);
            
        else
            dZ(iElec)=nan;
            Z(iElec,:)=nan(1,2);
        end
        
    end
    
    %take the mean of all times this channel is used
    Zave=mean(Z,2);
    
    %date stamp for this z check
    datestart=HDR.T0;
    datestart(6)=datestart(6)+Data_start;
    datestart=datenum(datestart);
    
    %% Classify channels
    %Each impedance is either good bad or ok
    bad_idx = find(Zave > Z_Max);
    ok_idx=find(Zave > Z_Rec & Zave < Z_Max);
    good_idx =find(Zave < Z_Rec);
    
    badchn=nan(size(Zave));
    goodchn=badchn;
    okchn=badchn;
    
    okchn(ok_idx)=Zave(ok_idx);
    goodchn(good_idx)=Zave(good_idx);
    badchn(bad_idx)=Zave(bad_idx);
    
    numbad=length(bad_idx);
    numok=length(ok_idx);
    
    %output to user
    disp('------------------------------------');
    fprintf('Found ');
    if numbad
        
        fprintf(2,'%d bad electrodes ',numbad);
        fprintf('and ');
    end
    
    fprintf('%d warning electrodes\n',numok);
    
    if numbad
        fprintf(2,'BAD ELECS : ');
        fprintf(2,'%d,',bad_idx);
        fprintf(2,'\n');
    end
    
    if numok
        fprintf('Warning elecs : ');
        fprintf('%d,',ok_idx);
        fprintf('\n');
    end
    
    
    %% plot results
    
    if PlotFlag
        
        figure
        title(sprintf('Z Check %d, in %s @ %s\n%d Hz and %d uA',iZ,fname,datestr(datestart),round(Fc{1}),Amp),'interpreter','none');
        
        hold all
        %plot each set
        bar(goodchn,'Facecolor',[0 0.5 0]);
        bar(okchn,'Facecolor','y');
        bar(badchn,'Facecolor',[1 0 0]);
        
        %make indication lines for recomended and max impedances
        line([0 Nelec+1],[Z_Max Z_Max],'color','r','linewidth',5)
        line([0 Nelec+1],[Z_Rec Z_Rec],'color','y','linewidth',5)
        text(1,Z_Rec,'FUZZY LOGIC OK','BackgroundColor',[1 1 1],'color',[0 0 0])
        text(1,Z_Max,'MAX Z','color','r','BackgroundColor',[1 1 1])
        hold off
        set(gca,'Xtick',[1:Nelec])
        xlim([0,Nelec+1])
        
        xlabel('Electrode');
        ylabel('~Impedance Ohm');
        
        
        %make plot wider
        pos=get(gcf,'Position');
        w=pos(1)-pos(3);
        set(gcf,'pos',[pos(1)-w pos(2) pos(3)+w pos(4)]);
        drawnow
    end
    
    %% Add to Structure
    
    %strucutre for a single z check
    Zout.Zi=Z;
    Zout.Z=Zave;
    Zout.dZ=dZ;
    Zout.ExpSetup=ExpSetup;
    Zout.info.Filt=Filt;
    Zout.info.Fc=Fc;
    Zout.info.FilterTrim=FilterTrim;
    Zout.info.bdf_filename=HDR.FILE.Name;
    Zout.TimeNum=datestart;
    Zout.TimeVec=datevec(datestart);
    
    %structure holding all
    Zall(iZ)=Zout;
end
%% save data
%find if the Zcheck directory exists by going up one directory and seeing
%if Zcheck exists.

Zcheckpathstr=HDR.FILE.Path;
bdfname=HDR.FILE.Name;

if Z_checks_num > 1
    
    Zfilename=[bdfname, '-Zcheck' num2str(Z_checks_num) '.mat'];
else
    Zfilename=[bdfname, '-Zcheck.mat'];
    
end

%save
save(fullfile(Zcheckpathstr,Zfilename),'Zall');
disp('----------------------');
end

