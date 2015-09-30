function [ BV,PhaseAngle ] = ScouseTom_ProcessBV( HDR,TT,ExpSetup,varargin )
%SCOUSETOM_ Summary of this function goes here
%   Detailed explanation goes here

% ONLY DOES ONE INEJCTION AT THE MOMENT

%% Do error checking of inputs HERE

%is there any data, do TT and ExpSetup match? Is the freq order ok?
%Expsetup and freq order should match

%allow for passing filter info into function


%% Defaults
BW=50; %bandwidth of bandpass filter in demod
progressbarlength=30; % number of - used in progressbar


%% See what type of system we are using


%the maximum voltage is different for each system, this *should* be in the
%HDR structure somewhere, but I dont know where it is

switch HDR.TYPE
    case 'BDF' % biosemi file
        
        MaxV=0.5e6; %500mV range on BioSemi
        
        
    case 'EEG'
end


%% get some variables from inputs Structures

Prot=ExpSetup.Protocol;
N_prt=size(Prot,1);
N_elec=ExpSetup.Elec_num;
N_freq=ExpSetup.Info.FreqNum;

Fs=HDR.SampleRate;
eegfname=HDR.FILE.Name;
eegfpath=HDR.FILE.Path;

%% Single or Multifrequency Mode

if N_freq == 1
    SingleFreqMode =1;
else
    SingleFreqMode=0;
end

%if in multifrequency mode then we need the freq order found in varargin

if ~SingleFreqMode
    FreqOrder=varargin{1};
end

%% calculate the keep and rem idx
% need this beore loading data to know which channels to estimate contact
% impedance on

[prt_full,keep_idx,rem_idx,Elec_inj]=ScouseTom_data_findprt(ExpSetup.Protocol,N_elec);

%scale factor - impedance conversion
ZSF=1/(ExpSetup.Amp); %keep this in uA as voltages are in uV

%% create matfile object for saving data

%matfile object can handle big files, and stops everything being stored in
%memory

% some basic info about the data set - HDR is too big to save each time
info.eegfname=eegfname;
info.TimeNum=datenum(HDR.T0);
info.TimeVec=HDR.T0;

%create matfile object in same place as data
bigmat=matfile(fullfile(eegfpath,[eegfname '-BV.mat']),'Writable',true);
% same tirgger info
bigmat.TT=TT;
%save the system settings
bigmat.ExpSetup=ExpSetup;
%save info and the protocol indexes
bigmat.info=info;
bigmat.keep_idx=keep_idx;
bigmat.rem_idx=rem_idx;
bigmat.prt_full=prt_full;

%store frewq order if it exists
if ~SingleFreqMode
    bigmat.FreqOrder=FreqOrder;
end


%% Calculate Data Length

%here we are cheating and are only taking first injection (thus covering
%99.99% of use cases!

%sread needs integer seconds
Data_start=floor(TT.InjectionStarts(1)/Fs);
Data_end=ceil(TT.InjectionStops(1)/Fs);
Data_length=fix(Data_end-Data_start);
Data_start_sample=Data_start*Fs;

Data_max=HDR.NRec; % max number of seconds - length of file

Data_start_s=Data_start*Fs;
Data_end_s=Data_end*Fs;


%% Get data in chunks!

% this bit is pretty badly written but i think it makes sense! chunks were
% taken in time rather than number of events as perhaps the injections were
% for a long time which would again cause out of memory errors

%start timer
tstart=tic;

%chunksize in seconds - how much data to load at once
chunk_time=10*60; %10 minutes

%finished loading flag
finished =0;
%first time loading flag
first=1;

%datawindow - the first seconds to be loaded
datawindow_s=[0 chunk_time];

%The Index of the switching vector last  used
next_sw=0;

%find the first line of the protocol - MAKE THIS CHECK RATHER THAN DO IT
lastprt=0;

% disp('----------------------');

%it
iteration=0;


curInjSwitch=TT.InjectionSwitches{1};
if ~SingleFreqMode
    curFreqSwitch=TT.FreqChanges{1};
else
    curFreqSwitch=[];
end


while finished == 0
    %% sort out which switches are in this chunk
    
    iteration=iteration+1;
    
    %convert window into samples to samples
    datawindow=(datawindow_s*Fs)+Data_start_sample;
    
    if first ==1
        %for first run find the start position
        idx_f=find(curInjSwitch > datawindow(1),1);
        %             first =0;
    else
        idx_f=next_sw;
    end
    % find last COMPLETE switch in recording - complete needs two switches
    idx_l=find(curInjSwitch > datawindow(2),1)-2;
    
    if isempty(idx_l) ==1;
        %if no more switches in file then this is the last iteration
        finished=1;
        idx_l=length(curInjSwitch);
    end
    %store this variable for next loop
    next_sw=idx_l+1;
    
    %find repeat - separate first case as var does not exist
    if first ==1
        start_rep=1;
    else
        start_rep=start_rep+N_rep;
    end
    
    
    %% load data in memory
    
    disp(['Loading data between ', num2str(datawindow_s(1)), 's and ', num2str(datawindow_s(2)), 's']);
    
    V=sread(HDR,datawindow_s(2)-datawindow_s(1),datawindow_s(1)+Data_start);
    %     V(:,N_elec+1:end)=[]; %remove extra channels in case reference was loaded too
    
    %%  Check Injection starts with expected injection
    
    %for first time only - find the injecting electrodes - query if it isnt the
    %first line in the protocol. This is to allow for files where the biosemi
    %stopped recording and we carried on.
    
    if first ==1
        [ lastprt ] = ScouseTom_data_checkfirstinj( V(1:Fs*2,:),Fs,Prot,curInjSwitch,curFreqSwitch,idx_f,datawindow,SingleFreqMode );
    end
    
    %next protocol line to use is one on from last one
    nextprt=lastprt+1;
    
    %% Find Filter Settings for each frequency
    
    if first ==1
        disp('--------Finding Filter Settings---------');
        
        if SingleFreqMode
            %using first injection
            tmp=curInjSwitch(idx_f+1)-curInjSwitch(idx_f);
            fwind=curInjSwitch(idx_f)-datawindow(1):curInjSwitch(idx_f)-datawindow(1)+tmp;
            
            %find carrier frequency and get filter coefficients as well as
            %the amount of data to remove each segment
            
            [trim_demod,B,A,Fc]=ScouseTom_data_GetFilterTrim(V(fwind,Prot(nextprt,1)),Fs,BW,0 );
        else
            %for multifreq do each one at a time
            for iFreq=1:N_freq
                
                f_idx_start=iFreq*2;
                f_idx_stop=f_idx_start+1;
                
                
                %take only the samples within the current frequency
                %injection
                fwind=(curFreqSwitch(f_idx_start):curFreqSwitch(f_idx_stop))-curInjSwitch(nextprt);
                
                
                [trim_demod{iFreq},B{iFreq},A{iFreq},Fc{iFreq}]...
                    =ScouseTom_data_GetFilterTrim( Vseg_demod(nextprt,fwind,Prot(nextprt,1),1),Fs,BW,0 );
            end
            
        end
        
        disp('--------Filter Settings Found------------');
        info.B=B;
        info.A=A;
        info.trim_demod=trim_demod;
        info.Fc=Fc;
        %save to .mat file
        bigmat.info=info;
    end
    
    %% Check Carrier frequecy is correct
    
    %do this here
    
    
    
    
    
    
    %% Demodulate Data - each channel at a time
    
    %demodulate data
    disp('Demodulating data');
    
    %put Voltages in matrix Sample x Channel x Frequency
    Vdemod=nan(size(V,1),size(V,2),N_freq);
    Pdemod=Vdemod; % this is matrix for Phase info
    
    %demodulate each freq in turn
    for iFreq=1:N_freq
        %demodulate entire channel at once
        for iElec=1:N_elec
            [Vdemod(:,iElec,iFreq),Pdemod(:,iElec,iFreq)] =ScouseTom_data_DemodHilbert(V(:,iElec),B(iFreq,:),A(iFreq,:));
        end
        
    end
    
    
    
    %% Segment data into chunks
    
    %multi freq will have to do this for each freq and put in cell array as
    %arrays are all different length
    
    
    
    %segment this data between the complete protocol lines
    [Vseg_demod,Pseg_demod, lastprt]=ScouseTom_data_Seg(Vdemod,Pdemod,curInjSwitch(idx_f:idx_l)-datawindow(1),0.0001,N_prt,N_elec,Fs,nextprt);
    
    %   If there are any left overs then stick them at the beginning
    if exist('Vsegleftover','var') ==1
        
        if size(Vseg_demod,2) == size(Vsegleftover,2)
            Vseg_demod(1:size(Vsegleftover,1),:,:,1)=Vsegleftover;
            Pseg_demod(1:size(Psegleftover,1),:,:,1)=Psegleftover;
        else if size(Vseg_demod,2) > size(Vsegleftover,2) %i cant remember why this is necessary - sometimes the chunks are different by a sample maybe?
                Vseg_demod(1:size(Vsegleftover,1),1:size(Vsegleftover,2),:,1)=Vsegleftover;
                Pseg_demod(1:size(Psegleftover,1),1:size(Psegleftover,2),:,1)=Psegleftover;
            else
                Vseg_demod(1:size(Vsegleftover,1),:,:,1)=Vsegleftover(:,1:size(Vseg_demod,2),:);
                Pseg_demod(1:size(Psegleftover,1),:,:,1)=Psegleftover(:,1:size(Pseg_demod,2),:);
            end
        end
        
        clear Vsegleftover
        clear Psegleftover
    end
    
    %take the incomplete repeat and trim matrix - only do this is there was
    %more than 1 repeat so Vseg is 4D
    if lastprt ~= N_prt && ndims(Vseg_demod) == 4
        Vsegleftover=Vseg_demod(1:lastprt,:,:,end);
        Vseg_demod(:,:,:,end)=[];
        Psegleftover=Pseg_demod(1:lastprt,:,:,end);
        Pseg_demod(:,:,:,end)=[];
    end
    
    %number of
    N_rep=size(Vseg_demod,4);
    N_sample=size(Vseg_demod,2);
    
    disp(['Number of complete repeats in chunk : ', num2str(N_rep)]);
    
    
    %% Calculate the standing BV for each channel
    disp('Getting Boundary Voltages');
    %get boundary voltages by taking mean
    [BV, STD]=ScouseTom_data_Seg2BV(Vseg_demod,trim_demod);
    %get phase angle by comparing to Injection channels
    [PhaseAngle,PhaseAngleSTD]=ScouseTom_data_PhaseEst(Pseg_demod,trim_demod,Prot);
    
    
    %% Calculate dZ and other Stimulation things
    
    %kirills code goes here
    
    
    
    %% Calculate Impedance on each injection electrode
    
    Z=nan(N_elec,N_rep);
    Zstd=Z;
    %get contact impedance values from BV
    for iElec=1:N_elec
        
        if all(isnan(Elec_inj(iElec,:)))
            Z(iElec,:)=nan;
            Zstd(iElec,:)=nan;
        else
            Z(iElec,:)=ZSF*nanmean(BV(Elec_inj(iElec,~isnan(Elec_inj(iElec,:))),:),1);
            Zstd(iElec,:)=ZSF*nanstd(BV(Elec_inj(iElec,~isnan(Elec_inj(iElec,:))),:),1);
        end
    end
    
    %% save Data
    
    %save them
    bigmat.BV(1:size(BV,1),start_rep:start_rep+N_rep-1)=BV;
    bigmat.STD(1:size(STD,1),start_rep:start_rep+N_rep-1)=STD;
    bigmat.Z(1:size(Z,1),start_rep:start_rep+N_rep-1)=Z;
    bigmat.Zstd(1:size(Zstd,1),start_rep:start_rep+N_rep-1)=Zstd;
    bigmat.PhaseAngle(1:size(PhaseAngle,1),start_rep:start_rep+N_rep-1)=PhaseAngle;
    bigmat.PhaseAngleSTD(1:size(PhaseAngle,1),start_rep:start_rep+N_rep-1)=PhaseAngleSTD;
    
    %% Output Progress
    
     timedone=toc(tstart);
    
    disp([num2str(timedone,'%.0f') 's:Finished processing between ', num2str(datawindow_s(1)), 's and ', num2str(datawindow_s(2)), 's']);
    
    %calculate the percentage complete
    percent_complete=100*(datawindow_s(2)/Data_max);
    
    if percent_complete >100
        percent_complete=100;
    end
    
    %display it visually
    
    decdone=floor(percent_complete/(100/progressbarlength));
    percentind=repmat('.',1,progressbarlength);
    percentind(1:decdone-1)=repmat('-',1,decdone-1);
    if decdone >0
        percentind(decdone)='>';
    end
    
    fprintf('%s %.1f %% complete\r',percentind,percent_complete);
    
    %% Calculate variables for next step
    
    %calculate new datablock
    
    if finished ==0
        
        %data window starts in second which includes next switch
        datawindow_s(1)=floor(curInjSwitch(next_sw)/Fs)-Data_start;
        %data window ends at next chunk
        datawindow_s(2)=datawindow_s(2)+chunk_time;
        
        if datawindow_s(2) >= Data_end
            datawindow_s(2) = Data_end;
        end
        
    end
    
    %this is no longer the first iteration!
    if first == 1
        first =0;
    end
    
    
end

%% All processing done!

disp('All processing finished! At f--ing last!');
teatime=toc(tstart);
fprintf('That took : %.1f seconds \r',teatime);

BV=bigmat.BV;




end

