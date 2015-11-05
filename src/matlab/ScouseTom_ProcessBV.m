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

%if in multifrequency mode then we need the freq order found in processing
%triggers

%% calculate the keep and rem idx
% need this beore loading data to know which channels to estimate contact
% impedance on

[prt_full,keep_idx,rem_idx,Elec_inj]=ScouseTom_data_findprt(ExpSetup.Protocol,N_elec);

%scale factor - impedance conversion
ZSF=1./(ExpSetup.Amp); %keep this in uA as voltages are in uV

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

%% Calculate Data Length

%here we are cheating and are only taking first injection (thus covering
%99.99% of use cases!

%sread needs integer seconds
Data_start=floor(TT.InjectionStarts(1)/Fs);
% Data_start=700;
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
chunk_time=5*60; %10 minutes
% % chunk_time=10;


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

%force first injection for now
curInjSwitch=[TT.InjectionSwitches{1};TT.InjectionStops(1)];
if ~SingleFreqMode
    curFreqOrder=TT.FreqOrder{1};
    curFreqStarts=TT.FreqStarts{1};
    curFreqStops=TT.FreqStops{1};
else
    curFreqOrder=[];
    curFreqStarts=[];
    curFreqStops=[];
end


while finished == 0
    %% sort out which switches are in this chunk
    
    if iteration == 5;
        disp('paasda');
    end
    
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
    % (also include stop switch here!
    idx_l=find([curInjSwitch]   > datawindow(2),1)-2;
    
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
        [ lastprt ] = ScouseTom_data_checkfirstinj( V(1:Fs*3,:),Fs,Prot,curInjSwitch,curFreqStarts,curFreqStops,idx_f,datawindow,SingleFreqMode );
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
            
            %make it consistent with multifreq bits, whic are all cells
            A={A};
            B={B};
            trim_demod={trim_demod};
        else
            %for multifreq do each one at a time
            for iFreq=1:N_freq
                
                freqidx=find(curFreqOrder(1,:) == iFreq);
                
                f_idx_start=curFreqStarts(1,freqidx);
                f_idx_stop=curFreqStops(1,freqidx);
                
                
                %take only the samples within the current frequency
                %injection
                fwind=(f_idx_start:f_idx_stop)-datawindow(1);
                
                [trim_demod{iFreq},B{iFreq},A{iFreq},Fc{iFreq}]...
                    =ScouseTom_data_GetFilterTrim( V(fwind,Prot(nextprt,1)),Fs,BW,0 );
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
            
            [Vdemod(:,iElec,iFreq),Pdemod(:,iElec,iFreq)] =ScouseTom_data_DemodHilbert(V(:,iElec),B{iFreq},A{iFreq});
        end
        
    end
    
    
    
    %% Segment data into chunks
    
    %Segment each dataset into chunks - taking freq is needed
    
    disp('Segmenting')
    
    %do this for each freq
    
    Vseg_demod=cell(N_freq,1);
    Pseg_demod=Vseg_demod;
    
    
    %% put data into cell array
    
    %im sorry for this data structure....
    for iFreq=1:N_freq
        
        %injection switches we need, wrt start of data NOT start of file
        Sw_seg=curInjSwitch(idx_f:idx_l)-datawindow(1);
        
        if SingleFreqMode
            Sw_seg=curInjSwitch(idx_f:idx_l)-datawindow(1);
            FreqStart_seg=[];
            FreqStop_seg=[];
        else
            
            
            %find the Starts and Stops related to this freq only
            FreqStart_seg_all=sort(curFreqStarts(curFreqOrder == iFreq)); % this is ALL of the frequency injections in file
            FreqStop_seg_all=sort(curFreqStops(curFreqOrder == iFreq));
            
            freq_idx_f=idx_f; %first index for freqs we can assume is ok
            
            %if injection was stopped early, the number of
            %complete injections could be different for each freq. I.e. ONly the ones
            %completed fully would give a FreqStop_Seg_all the maximum length idx_l-1, so we
            %have to find the maximum index separately, as
            
            if idx_l-1 > size(FreqStart_seg_all,1)
                
                
                freq_idx_l=size(FreqStart_seg_all,1);
            else
                freq_idx_l=idx_l-1;
            end
            
            
            %only want the ones within the data loaded
            FreqStart_seg=FreqStart_seg_all(freq_idx_f:freq_idx_l)-datawindow(1);
            FreqStop_seg=FreqStop_seg_all(freq_idx_f:freq_idx_l)-datawindow(1);
            
            Sw_seg=curInjSwitch(freq_idx_f:idx_l)-datawindow(1);
            
        end
        %segment this data between the complete protocol lines
        [Vseg_demod{iFreq},Pseg_demod{iFreq}, lastprt]=ScouseTom_data_Seg(Vdemod(:,:,iFreq),Pdemod(:,:,iFreq),Sw_seg,FreqStart_seg,FreqStop_seg,0.0001,N_prt,N_elec,Fs,nextprt);
        
    end
    %% Any left overs?
    
    %   If there are any left overs then stick them at the beginning
    if exist('Vsegleftover','var') ==1
        
        for iFreq=1:N_freq
            
            if size(Vseg_demod{iFreq},2) == size(Vsegleftover{iFreq},2)
                Vseg_demod{iFreq}(1:size(Vsegleftover{iFreq},1),:,:,1)=Vsegleftover{iFreq};
                Pseg_demod{iFreq}(1:size(Psegleftover{iFreq},1),:,:,1)=Psegleftover{iFreq};
            else if size(Vseg_demod{iFreq},2) > size(Vsegleftover{iFreq},2) %i cant remember why this is necessary - sometimes the chunks are different by a sample maybe?
                    Vseg_demod{iFreq}(1:size(Vsegleftover{iFreq},1),1:size(Vsegleftover{iFreq},2),:,1)=Vsegleftover{iFreq};
                    Pseg_demod{iFreq}(1:size(Psegleftover{iFreq},1),1:size(Psegleftover{iFreq},2),:,1)=Psegleftover{iFreq};
                else
                    Vseg_demod{iFreq}(1:size(Vsegleftover{iFreq},1),:,:,1)=Vsegleftover{iFreq}(:,1:size(Vseg_demod{iFreq},2),:);
                    Pseg_demod{iFreq}(1:size(Psegleftover{iFreq},1),:,:,1)=Psegleftover{iFreq}(:,1:size(Pseg_demod{iFreq},2),:);
                end
            end
            
        end
        
        
        clear Vsegleftover
        clear Psegleftover
    end
    
    %make leftover cell for each freq
    for iFreq=1:N_freq
        
        %take the incomplete repeat and trim matrix - only do this is there was
        %more than 1 repeat so Vseg is 4D
        if lastprt ~= N_prt && ndims(Vseg_demod{iFreq}) == 4
            Vsegleftover{iFreq}=Vseg_demod{iFreq}(1:lastprt,:,:,end);
            Vseg_demod{iFreq}(:,:,:,end)=[];
            Psegleftover{iFreq}=Pseg_demod{iFreq}(1:lastprt,:,:,end);
            Pseg_demod{iFreq}(:,:,:,end)=[];
        end
    end
    
    %% output to user
    
    %number of
    N_rep=size(Vseg_demod{1},4);
    
    disp(['Number of complete repeats in chunk : ', num2str(N_rep)]);
    
    
    %% Calculate the standing BV for each channel
    disp('Getting Boundary Voltages');
    
    %     %Vdemod is not a cell if only 1 freq
    %     if SingleFreqMode
    %
    %          [BV, STD]=ScouseTom_data_Seg2BV(Vseg_demod,trim_demod);
    %         %get phase angle by comparing to Injection channels
    %         [PhaseAngle,PhaseAngleSTD]=ScouseTom_data_PhaseEst(Pseg_demod,trim_demod,Prot);
    %
    %
    %
    %     else
    
    
    clear BV STD PhaseAngle PhaseAngleSTD
    for iFreq=1:N_freq
        
        %get boundary voltages by taking mean
        [BV(:,:,iFreq), STD(:,:,iFreq)]=ScouseTom_data_Seg2BV(Vseg_demod{iFreq},trim_demod{iFreq});
        %get phase angle by comparing to Injection channels
        [PhaseAngle(:,:,iFreq),PhaseAngleSTD(:,:,iFreq)]=ScouseTom_data_PhaseEst(Pseg_demod{iFreq},trim_demod{iFreq},Prot);
        
    end
    
    %     end
    %% Calculate dZ and other Stimulation things
    
    %kirills code goes here
    
    
    %% clear variables
    
    %     clear Vseg_demod Pseg_demod Vdemod Pdemod
    
    
    
    %% Calculate Impedance on each injection electrode
    
    Z=nan(N_elec,N_rep,N_freq);
    Zstd=Z;
    
    for iFreq=1:N_freq
        
        %get contact impedance values from BV
        for iElec=1:N_elec
            
            %set to nan is there were no injections on this elec
            if all(isnan(Elec_inj(iElec,:)))
                Z(iElec,:,iFreq)=nan;
                Zstd(iElec,:,iFreq)=nan;
            else %average all the voltages on the injection channel
                Z(iElec,:,iFreq)=ZSF(iFreq)*nanmean(BV(Elec_inj(iElec,~isnan(Elec_inj(iElec,:))),:,iFreq),1);
                Zstd(iElec,:,iFreq)=ZSF(iFreq)*nanstd(BV(Elec_inj(iElec,~isnan(Elec_inj(iElec,:))),:,iFreq),1);
            end
        end
        
    end
    
    %% save Data
    
    %save them to the matfile
    
    if SingleFreqMode %have to do this differently as the matfile doesnt like 3 indicies even if N_freq is 1
        
        bigmat.BV(1:size(BV,1),start_rep:start_rep+N_rep-1)=BV;
        bigmat.STD(1:size(STD,1),start_rep:start_rep+N_rep-1)=STD;
        bigmat.Z(1:size(Z,1),start_rep:start_rep+N_rep-1)=Z;
        bigmat.Zstd(1:size(Zstd,1),start_rep:start_rep+N_rep-1)=Zstd;
        bigmat.PhaseAngle(1:size(PhaseAngle,1),start_rep:start_rep+N_rep-1)=PhaseAngle;
        bigmat.PhaseAngleSTD(1:size(PhaseAngle,1),start_rep:start_rep+N_rep-1)=PhaseAngleSTD;
        
        
    else
        bigmat.BV(1:size(BV,1),start_rep:start_rep+N_rep-1,1:N_freq)=BV;
        bigmat.STD(1:size(STD,1),start_rep:start_rep+N_rep-1,1:N_freq)=STD;
        bigmat.Z(1:size(Z,1),start_rep:start_rep+N_rep-1,1:N_freq)=Z;
        bigmat.Zstd(1:size(Zstd,1),start_rep:start_rep+N_rep-1,1:N_freq)=Zstd;
        bigmat.PhaseAngle(1:size(PhaseAngle,1),start_rep:start_rep+N_rep-1,1:N_freq)=PhaseAngle;
        bigmat.PhaseAngleSTD(1:size(PhaseAngle,1),start_rep:start_rep+N_rep-1,1:N_freq)=PhaseAngleSTD;
    end
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

