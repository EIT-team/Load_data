function [ BV ] = ScouseTom_ProcessBV( HDR,TT,ExpSetup )
%SCOUSETOM_ Summary of this function goes here
%   Detailed explanation goes here

% ONLY DOES ONE INEJCTION AT THE MOMENT

%% Do error checking of inputs HERE

%is there any data, do TT and ExpSetup match? 



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

Fs=HDR.SampleRate;

eegfname=HDR.FILE.Name;
eegfpath=HDR.FILE.Path;


BW=50; %bandwidth of bandpass filter in demod


%% calculate the keep and rem idx
% need this beore loading data to know which channels to estimate contact
% impedance on

injs=ExpSetup.Protocol; %injection pairs
chn=N_elec; %number of channels
vp=(1:chn)';% positive voltage channel 1  to Number of electrodes
vm=ones(size(vp))*(chn+1); %always against a ground electrode

%make the entire protocol each line is INJ+ INJ- MEAS+ MEAS-
prt_mat=[];
for iii=1:size(injs,1)
    temp=[repmat(injs(iii,:),chn,1) vp vm];
    prt_mat=[prt_mat ; temp];
end

%find remove index, any protocol lines including the injetion channels are
%"bad"
prt_full=prt_mat;
prt=prt_full;
rem_idx=[];
for iPrt = 1:size(prt,1)
    if any(ismember(prt_full(iPrt,1:2),prt(iPrt,3:4))) ==1
        rem_idx=[rem_idx,iPrt];
    end
end
%keep index is anything that we *dont* remove
keep_idx=setdiff(1:length(prt_full),rem_idx);

%% get injection channels for use in contact impedance calculations

%loop through protocol - find which lines in the BV the injection
%electrodes belong to then add then to an array for each electrode.

%I cant remember why I do this separately to the stuff about injections
%above...

%Electrode Injections
Elec_inj=nan(N_elec,N_prt);

for iPrt = 1:N_prt
    Prt_cur=Prot(iPrt,:);
    start_idx=((iPrt-1)*N_elec);
    BV_chn=start_idx+Prt_cur;
    Elec_inj(Prt_cur,iPrt)=BV_chn;
end

Elec_inj=sort(Elec_inj,2);

%clear up matrix
Elec_inj(:,all(isnan(Elec_inj),1))=[];

%scale factor - impedance conversion
ZSF=1/((1e6)*ExpSetup.Amp);

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

disp('----------------------');

%it
iteration=0;


curInjSwitch=TT.InjectionSwitches{1};


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
    
    % find last COMPLETE switch in recording
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
        
        %take either the first injection or the first second
        
        tmp=curInjSwitch(idx_f+1)-curInjSwitch(idx_f);
        if tmp > Fs
            tmp=Fs;
        end
        
        tmpidx=curInjSwitch-datawindow(1):curInjSwitch-datawindow(1)+tmp;
        
        %estimate the injection pairs from the two largest RMS values
        [InjPairs, estimatebadness]=ScouseTom_data_EstInjPair(V(tmpidx,:));
        
        %get the injection pair from the protocol
        ProtPairs=Prot(1,:)';
        
        %if the estimation is OK then calculate automatically
        
        if estimatebadness == 0
            
            %if the injection channels match then crack on
            if all(sort(InjPairs)==sort(ProtPairs)) ==1
                %disp('Data starts with first protocol line');
                lastprt=0; %this is already set above but being didactic
            else
                disp('----------------------');
                disp('Data DOES NOT start with first protocol line');
                
                %find matching protocl line
                start_poss=find(all([InjPairs(1)==Prot(:,1) InjPairs(2)==Prot(:,2)],2));
                disp(['Starting injection pair was found to be : ', num2str(start_poss)])
                disp('Data processing carrying on now...');
                
                lastprt=start_poss-1;
            end
            
        else
            disp('----------------------');
            %if it is still ambiguous - ask the user what to do
            msgbox('Starting injection pair is ambiguous! Please check the graph and enter manually','Uh Oh!');
            
            %plot the voltages for this swithc
            figure;
            plot(V(tmpidx,:));
            title('starting injection data - ambiguous injection');
            
            %ask them to input which protocol line this
            start_poss=input('Please enter the protocol line or leave empty to use best guess:');
            
            %is its empty just use best guess
            if isempty(start_poss)
                disp('FINE! I will just use the possibly wrong guess then shall I?');
                start_poss=find(all([InjPairs(1)==Prot(:,1) InjPairs(2)==Prot(:,2)],2));
                
                disp(['Starting injection pair was found to be : ', num2str(start_poss)'])
                disp('Data processing carrying on now...');
            end
            lastprt=start_poss-1;
            disp('----------------------');
        end
        
        
    end
    
    %% Segment data into each injection
    
%     disp('Segmenting data');
    
    %next protocol line to use is one on from last one
    nextprt=lastprt+1;
    
    %segment this data between the complete protocol lines starting
    %from correct protocol line
    [Vseg, lastprt]=ScouseTom_data_Seg(V,curInjSwitch(idx_f:idx_l)-datawindow(1),0.0001,N_prt,N_elec,Fs,nextprt);
    
    % for the first time only - determine carrier frequency, filter coeffs,
    % samples to trim in demodulation
    
    if first ==1
        disp('----------------------');
        %using first injection, find the best filter coefficients
        [trim_demod,B,A,Fc]=ScouseTom_data_GetFilterTrim( Vseg(nextprt,:,Prot(nextprt,1),1),Fs,BW,0 );
        disp('----------------------');
        info.B=B;
        info.A=A;
        info.trim_demod=trim_demod;
        info.Fc=Fc;
        %save to .mat file
        bigmat.info=info;
    end
    
    %   If there are any left overs then stick them at the beginning
    
    
    if exist('Vsegleftover','var') ==1
        
        if size(Vseg,2) == size(Vsegleftover,2)
            Vseg(1:size(Vsegleftover,1),:,:,1)=Vsegleftover;
        else if size(Vseg,2) > size(Vsegleftover,2)
                Vseg(1:size(Vsegleftover,1),1:size(Vsegleftover,2),:,1)=Vsegleftover;
            else
                Vseg(1:size(Vsegleftover,1),:,:,1)=Vsegleftover(:,1:size(Vseg,2),:);
            end
        end
        
        clear Vsegleftover
    end
    
    
    %take the incomplete repeat and trim matrix - only do this is there was
    %more than 1 repeat so Vseg is 4D
    if lastprt ~= N_prt && ndims(Vseg) == 4
        Vsegleftover=Vseg(1:lastprt,:,:,end);
        Vseg(:,:,:,end)=[];
    end
    
    %number of
    N_rep=size(Vseg,4);
    N_sample=size(Vseg,2);
    
    disp(['Number of complete repeats in chunk : ', num2str(N_rep)]);
    
    %% Demodulate and calculate BV from average in each segment
    
    
    %demodulate data
    disp('Demodulating data');
    Vseg_demod=ScouseTom_data_DemodSeg(Vseg,Fs,N_prt,N_elec,N_rep,B,A);
    
    disp('Getting Boundary Voltages');
    
    %get boundary voltages by taking mean
    [BV, STD]=ScouseTom_data_Seg2BV(Vseg_demod,trim_demod);
    %     clear V_seg_demod
    
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
    
    %% Output Progress
    disp(['Finished processing data between ', num2str(datawindow_s(1)), 's and ', num2str(datawindow_s(2)), 's']);
    
    %calculate the percentage complete
    percent_complete=100*(datawindow_s(2)/Data_max);
    
    if percent_complete >100
        percent_complete=100;
    end
    
    %display it visually
    winlength=30;
    decdone=floor(percent_complete/(100/winlength));
    percentind=repmat('.',1,winlength);
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

