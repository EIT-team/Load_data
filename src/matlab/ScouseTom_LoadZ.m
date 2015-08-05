function [ Zall ] = ScouseTom_LoadZ( HDR,TT,ExpSetup,PlotFlag )
%SCOUSETOM_LOADZ Calculates the contact impedance from a two electrode "contact check"
%   Detailed explanation goes here

%% check inputs are ok here

if exist('PlotFlag','var') ==0 
    PlotFlag =1;
end




%% See what type of system we are using


%the maximum voltage is different for each system, this *should* be in the
%HDR structure somewhere, but I dont know where it is

switch HDR.TYPE
    case 'BDF' % biosemi file
        
        MaxV=0.5e6; %500mV range on BioSemi
        
        
    case 'EEG'
end

%% Get variables

Fs=HDR.SampleRate;

N_elec=ExpSetup.Elec_num; %number of electrodes
N_prt=N_elec; %this is the same as N_elec as the protocol lines is equal to number of electrodes for contact checks



%% Amplitude and freq is FIXED AT THE MOMENT
Amp=141; % 141 uA
Frq=1000;



%% Defaults

%bandwidth of filter
BW=50;
RecZ=1000; %recommended contact impedance


%% run each impedance check separately

% run each
Z_checks_num=length(TT.InjectionStarts);

for iZ = 1:Z_checks_num
    
    
    %% read portion of interest
    
    %sread needs integer seconds
    Data_start=floor(TT.InjectionStarts(iZ)/Fs);
    Data_end=ceil(TT.InjectionStops(iZ)/Fs);
    Data_length=fix(Data_end-Data_start);
    
    %convert to number of samples
    Data_start_s=Data_start*Fs;
    Data_end_s=Data_end*Fs;
    
    %load the data
    disp(['Reading relevant data for contact check :', num2str(iZ)]);
    
    %sread MUST take INTEGER seconds as inputs, hence rounding data to
    %nearest second as above
    V=sread(HDR,Data_length,Data_start);
    
    % data is now in big long streams - segment into separate lines of the
    % protocol - output is matrix of voltages size PRT x Sample x CHN x Repeat
    Vseg=ScouseTom_data_Seg(V(:,1:N_elec),TT.InjectionSwitches{iZ}-Data_start_s,0.0001,N_prt,N_elec,Fs);
    
    %determine the Carrier Frequency, Filter coeffs, and amount of signal
    %to trim from the electrode on the *second* injection
    [trim_demod,B,A,Fc]=ScouseTom_data_GetFilterTrim( Vseg(2,:,2),Fs,BW,0 );
    
    
    %demodulate each segment in turn using hilber transfrom
    fprintf('Demodulating...');
    Vsegdemod=ScouseTom_data_DemodSeg(Vseg,Fs,N_prt,N_elec,1,B,A);
    fprintf('Done!\n');
    
    %get boundary voltages - average each segment - ignoring the messed up
    %edges
    BV=ScouseTom_data_Seg2BV(Vsegdemod,trim_demod,'please dont reshape thank youu please');
    
    %only interested in the injection voltages
    
    Injchn=[1:N_elec; circshift(1:N_elec,[1 -1])]';
    InjchnV=nan(size(Injchn));
    Injchndiff=nan(N_elec,1);
    
    
    %get the injection voltages from each pair as well as the "legacy"
    %difference measurements
    for iElec=1:N_elec;
        
        InjchnV(iElec,:)=BV(iElec,Injchn(iElec,:));
        Injchndiff(iElec)=diff(InjchnV(iElec,:));
    end
    
    %convert into impedances using peak amplitude values and convert uV to V
    
    %use peak values as peak values calculated during demod
    
    %scale factor - impedance conversion
    ZSF=1/(Amp); %VOLTAGE GIVEN IN uV SO AMP in uA too
    
    MaxZ=MaxV*ZSF;
    
    %calc impedances and align each repeat for electrode
    Z=InjchnV*ZSF;
    Z=[Z(:,1), circshift(Z(:,2),[1 1])];
    %calc the "legacy" UCH 2 channel one too
    dZ=Injchndiff*ZSF;
    
    
    %take average of the two readings for each electrode
    Zave=nanmean(Z,2);
    
    %date stamp for this z check
    datestart=HDR.T0;
    datestart(6)=datestart(6)+Data_start;
    datestart=datenum(datestart);
    
    %% plot results
    
    if PlotFlag
    
    

       
    figure
    
    title([datestr(datestart), ': Contact Z @ ', num2str(Fc), ' Hz, and ', num2str(Amp), ' A'])
    hold all
    
    
    bad_idx = find(Zave > MaxZ);
    ok_idx=find(Zave > RecZ & Zave < MaxZ);
    good_idx =find(Zave < RecZ);
    
   badchn=nan(size(Zave));
   goodchn=badchn;
   okchn=badchn;
   
   okchn(ok_idx)=Zave(ok_idx);
   goodchn(good_idx)=Zave(good_idx);
   badchn(bad_idx)=Zave(bad_idx);
   
    
   
    
    bar(goodchn,'Facecolor',[0 0.5 0]);
    bar(okchn,'Facecolor','y');
    bar(badchn,'Facecolor',[1 0 0]);
    
   
    line([0 N_elec+1],[MaxZ MaxZ],'color','r','linewidth',5)
    line([0 N_elec+1],[RecZ RecZ],'color','y','linewidth',5)
    text(1,RecZ,'FUZZY LOGIC OK','BackgroundColor',[1 1 1],'color',[0 0 0])
    text(1,MaxZ,'MAX Z','color','r','BackgroundColor',[1 1 1])
    hold off
    set(gca,'Xtick',[1:N_elec])
    xlim([0,N_elec+1])
    
    xlabel('Electrode');
    ylabel('~Impedance Ohm');
    
    
    end
    
    Zout.Zi=Z;
    Zout.Z=Zave;
    Zout.dZ=dZ;
    Zout.ExpSetup=ExpSetup;
    Zout.info.B=B;
    Zout.info.A=A;
    Zout.info.Fc=Fc;
    Zout.info.trim_demod=trim_demod;
    Zout.info.bdf_filename=HDR.FILE.Name;
    Zout.TimeNum=datestart;
    Zout.TimeVec=datevec(datestart);
    
    Zall(iZ)=Zout;
    
    
    
    
    
    
end
%% save data
    %find if the Zcheck directory exists by going up one directory and seeing
    %if Zcheck exists.
    
Zcheckpathstr=HDR.FILE.Path;
bdfname=HDR.FILE.Name;

        if Z_checks_num > 1
            
            Zfilename=[bdfname, '-Zcheck' num2str(iZ) '.mat'];
        else
            Zfilename=[bdfname, '-Zcheck.mat'];

        end
    
    
    %save
    save(fullfile(Zcheckpathstr,Zfilename),'Zall');
    disp('----------------------');
    


end

