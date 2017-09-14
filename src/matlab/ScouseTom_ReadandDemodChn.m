function [ VmagOut,PhaseOut,VmagOutSTD,PhaseOutSTD,Vmag ] = ScouseTom_ReadandDemodChn( HDRin,Filt,Trim_demod,InjectionWindows,Protocol,StartInj )
% [VmagOut,PhaseOut,VmagOutSTD,PhaseOutSTD,Vmag ] = ScouseTom_ReadandDemodChn( HDRin,Filt,Trim_demod,InjectionWindows,Protocol,StartInj )

%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


%% Max memory usage
MaxMemoryUsage=8e9; %maximum memory usage for larger variables stored during demodulation

%% Read HDR and get file size
HDR=HDRin; %backup HDR
switch HDR.TYPE
    case 'BDF' % biosemi file
        %biosemi is stored as 1sec long records. so number of seconds is
        %Nrec
        Nsec=HDR.NRec;
        %biosemi file size is stored in bytes
        Fsize=HDR.FILE.size;
    case 'BrainVision'
        Nsec=ceil(HDR.SPR/HDR.SampleRate);
        %actichamp HDR points to header .vhdr file, so filesize is wrong
        eegfile=dir(fullfile(HDR.FILE.Path,[HDR.FILE.Name '.eeg']));
        Fsize=eegfile.bytes;
    otherwise
        error('Bad HDR');
end
Nchn=size(HDR.InChanSelect,1);
Fs=HDR.SampleRate;

%% Check inputs

Nfreq=size(InjectionWindows,2);
Nprt=size(Protocol,1);

% Assume we start at the beginning if not specified
if exist('StartInj','var') ==0
    StartInj=ones(Nfreq,1);
    fprintf(2,'WARNING. NO STARTINJ SPECIFIED. ASSUMING 1. \n');
end
% check if we have start injection pair for all freqs
if length(StartInj) ~= Nfreq
    error('Wrong number of start injections');
end

%% Find start and finish times

% Start and stop times are given in TT as samples, but the sread function needs
% secons, and these need to be integers.

firstinj=zeros(Nfreq,1);
lastinj=zeros(Nfreq,1);
%find earliest switch
for iFreq=1:Nfreq
    firstinj(iFreq)=min(InjectionWindows{iFreq}(1,:));
    lastinj(iFreq)=max(InjectionWindows{iFreq}(end,:));
end

[FirstSample]=min(firstinj);
[LastSample]=max(lastinj);

%find the start time in seconds - subtract 1 to give more samples to reduce
%filtering artefacts
StartSec=floor(FirstSample/Fs)-1;
if StartSec < 0
    StartSec =0; %make sure it is >=0
end

%find the stop time in seconds - add 1 to give more samples to reduce
%filtering artefacts
StopSec=ceil(LastSample/Fs)+1;
if StopSec > Nsec
    StopSec=StopSec-1; %remove added second if we go too far
end

Start_Sample=StartSec*Fs;

if StopSec > Nsec
    fprintf(2,'Stop sec is too long!\n');
end

%% preallocate

Vmag=cell(Nfreq,1);

for iFreq=1:Nfreq
    Vmagtmp=nan(size(InjectionWindows{iFreq},1),Nchn);
    Vmag{iFreq}=Vmagtmp;
end
PhaseRaw=Vmag;
VmagSTD=Vmag;
PhaseRawSTD=Vmag;
Phase=Vmag;

%% Based on file size, determine how to load data
% To limit the number of times we read the data, we want to demodulate as
% many channels as possible, but we run out of memory. So this estimates
% how much RAM is needed per channel, and calculates the maximum number we
% can fit.
% This is hacky, sorry!

%file size in bytes is ~double when stored in matlab as double
ChnSize=(Fsize*2.4)/Nchn;
%this is because we have to store the demod V and Phase for each channel -
%add more because the fft in hilbert takes LOADS of ram
ChnSize=ChnSize*4;

%maximum number of channels this length which could be put into memory
MaxChnNum=floor(MaxMemoryUsage/ChnSize);

%flag to see whether we need to load channel by segment or not. If a
%channel is so big - like Nirs 7 hour recordings - then even a single
%channel cannot be processed at once. So we have to treat these cases
%differently
LoadSegmentsFlag = 0;

%force us to use at least one channel at once, this could still mean you
%run out of memory later as whole channel is still loaded!
if MaxChnNum <= 1
    MaxChnNum =1; %need this to be 1 so we actually load a channel
    LoadSegmentsFlag = 1; % flag up we need to load it differently
    NumSeg = 4; % How many segments to split the channels into. might need to calc this separately later...
    fprintf('Loading Single Channel. Processing in %d segments to avoid memory issues.\n', NumSeg);
end
% how many blocks do we have to split channels into
BlocksNum=ceil(Nchn/MaxChnNum);
%ending channel for each block
Blocks=((1:BlocksNum)).*MaxChnNum;

%create the array of start and stop channel pairs.

if MaxChnNum == 1
    Blocks=[Blocks; Blocks]'; % [1 1; 2 2;...]
else
    %starting channel is ending minus maxchnnum then adjusted for 1 indexing
    Blocks=[Blocks-(MaxChnNum-1); Blocks]';
    %to prevent loading channels that dont exist
    Blocks(Blocks > Nchn)=Nchn;
end
%% Read and Demod each channel

fprintf('Processing %d Channels',Nchn);
if Nfreq > 1
    fprintf(' at %d Frequencies',Nfreq);
end
fprintf('\n');

tstart=tic;

% read and demodulate the data, taking each block of channels at a time
for iBlk=1:BlocksNum
    
    fprintf('Process Chn %d to %d. Freq: ',Blocks(iBlk,1),Blocks(iBlk,2));
    
    curChn=Blocks(iBlk,1):Blocks(iBlk,2);
    curChnNum=size(curChn,2);
    
    %set variables in HDR for single channel only
    HDR.InChanSelect=curChn;
    HDR.Calib=HDRin.Calib(1:curChnNum+1,1:curChnNum);
    
    V=sread(HDR,StopSec-StartSec,StartSec); %read whole channel
    
    %% Demodulate this data
    %demodulate for each frequency
    for iFreq=1:Nfreq
        %display which freq we are doing
        if iFreq < Nfreq
            fprintf('%d,',iFreq);
        else
            fprintf('%d',iFreq);
        end
        
        % If we have determined that a single channel is too big to
        % process at once, then do the filtering and demodulation in
        % segments
        if LoadSegmentsFlag
            %split single channel into segments
            Vdata_demod = nan(size(V));
            Pdata_demod = nan(size(V));
            
            % where we are going to split channel
            Seg_Samples = [ ceil(1:length(V)/NumSeg:length(V)), length(V)];
            
            % we need to add buffers to avoid edge effects from the
            % filtering
            
            BufferSamples = 10*Fs;
            Seg_Start = [Seg_Samples(1), Seg_Samples(2:end-1) - BufferSamples];
            Seg_End= [Seg_Samples(2:end-1) + BufferSamples, length(V)];
            
            %filter and demodulate per segment
            for iSeg=1:NumSeg
                
                %get filtered hibs from this data segment
                [Vtmp,Ptmp] = ScouseTom_data_DemodHilbert( V(Seg_Start(iSeg):Seg_End(iSeg)),Filt{iFreq});
                
                %find the bit within this segment we want
                cur_idx_start = max([Seg_Samples(iSeg) - Seg_Start(iSeg)+1, 1]);
                cur_idx_stop = min([length(Vtmp),cur_idx_start+(Seg_Samples(iSeg+1)-Seg_Samples(iSeg))]);
                %stick it in the full dataset
                Vdata_demod(Seg_Samples(iSeg):Seg_Samples(iSeg+1)) = Vtmp(cur_idx_start:cur_idx_stop);
                Pdata_demod(Seg_Samples(iSeg):Seg_Samples(iSeg+1))  = Ptmp(cur_idx_start:cur_idx_stop);
            end
                        
        else
            
            % filter and demodulate channels all at once
            [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,Filt{iFreq});
            
        end
        %process each injection window, adjusting for new start time
        [Vmag{iFreq}(:,curChn),PhaseRaw{iFreq}(:,curChn),VmagSTD{iFreq}(:,curChn),PhaseRawSTD{iFreq}(:,curChn)]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod{iFreq},InjectionWindows{iFreq}-Start_Sample);
    end
    
    t_el=toc(tstart);
    fprintf('. t=%.1f s\n',t_el);
    
end

teatime=toc(tstart);
fprintf('ALL DONE! That took : %.1f seconds \r',teatime);

%% Calculate Phase

%phase is absolute (or rather with respect to sin/cos starting at sample 0
%of whole dataset. So we need to get the phase relative to the injection
%electrodes
for iFreq=1:Nfreq
    [Phase{iFreq}]=ScouseTom_data_PhaseEst(PhaseRaw{iFreq},Protocol,StartInj);
end

%% Pad to complete number of repeats

for iFreq=1:Nfreq
    
    %number of injections found
    Ninj=size(Vmag{iFreq},1);
    
    %find number of repeats - rounding up
    Nrep=ceil((StartInj(iFreq)+Ninj)/Nprt);
    
    %preallocate
    VmagOutTmp=nan(Nchn,Nrep*Nprt);
    PhaseOutTmp=VmagOutTmp;
    VmagOutSTDTmp=VmagOutTmp;
    PhaseOutSTDTmp=VmagOutTmp;
    
    %put data we have in the correct place
    VmagOutTmp(:,StartInj(iFreq):(StartInj(iFreq)-1)+Ninj)=Vmag{iFreq}';
    PhaseOutTmp(:,StartInj(iFreq):(StartInj(iFreq)-1)+Ninj)=Phase{iFreq}';
    VmagOutSTDTmp(:,StartInj(iFreq):(StartInj(iFreq)-1)+Ninj)=VmagSTD{iFreq}';
    PhaseOutSTDTmp(:,StartInj(iFreq):(StartInj(iFreq)-1)+Ninj)=PhaseRawSTD{iFreq}';
    
    %% reshape into correct format
    VmagOut{iFreq}=reshape(VmagOutTmp,Nchn*Nprt,Nrep);
    PhaseOut{iFreq}=reshape(PhaseOutTmp,Nchn*Nprt,Nrep);
    VmagOutSTD{iFreq}=reshape(VmagOutSTDTmp,Nchn*Nprt,Nrep);
    PhaseOutSTD{iFreq}=reshape(PhaseOutSTDTmp,Nchn*Nprt,Nrep);
    %% remove empty repeats
    %easier to it this way that to ensure it doesnt happen higher up
    
    VmagOut{iFreq}(:,all(isnan(VmagOut{iFreq})))=[];
    PhaseOut{iFreq}(:,all(isnan(PhaseOut{iFreq})))=[];
    VmagOutSTD{iFreq}(:,all(isnan(VmagOutSTD{iFreq})))=[];
    PhaseOutSTD{iFreq}(:,all(isnan(PhaseOutSTD{iFreq})))=[];
    
end

end

