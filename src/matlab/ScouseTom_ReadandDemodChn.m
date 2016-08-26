function [ VmagOut,PhaseOut,VmagOutSTD,PhaseOutSTD,Vmag ] = ScouseTom_ReadandDemodChn( HDRin,Filt,Trim_demod,InjectionWindows,Protocol,StartInj )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Fiddling with inputs
HDR=HDRin;
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

if exist('StartInj','var') ==0
    StartInj=ones(Nfreq,1);
    fprintf(2,'WARNING. NO STARTINJ SPECIFIED. ASSUMING 1. \n');
end

if length(StartInj) ~= Nfreq
    error('Wrong number of start injections');
end

%% Find start and finish times

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
MaxMemoryUsage=8e9; %maximum memory usage for larger variables stored during demodulation
% MaxMemoryUsage=.000001e9; %maximum memory usage for V variable in bytes

ChnSize=(Fsize*2.4)/Nchn; %file size in bytes is ~double when stored in matlab as double

ChnSize=ChnSize*4; %this is because we have to store the demod V and Phase for each channel - add more because the fft in hilbert takes LOADS of ram

MaxChnNum=floor(MaxMemoryUsage/ChnSize); %maximum number of channels this length which could be put into memory

%force us to use at least one channel at once, this could still mean you
%run out of memory later!
if MaxChnNum <= 1
    MaxChnNum =1; %need this to be 1 so we actually load a channel
    LoadSegmentsFlag = 1;
    
end
BlocksNum=ceil(Nchn/MaxChnNum); %how many blocks do we have to split channels into

Blocks=((1:BlocksNum)).*MaxChnNum; %ending channel for each block

if MaxChnNum == 1
    
    
     Blocks=[Blocks; Blocks]';
     
     
else
    
    
    
    Blocks=[Blocks-(MaxChnNum-1); Blocks]'; %starting channel is ending minus maxchnnum then adjusted for 1 indexing
    
    Blocks(Blocks > Nchn)=Nchn; %to prevent loading channels that dont exist
end



%% Read and Demod each channel

fprintf('Processing %d Channels',Nchn);
if Nfreq > 1
    fprintf(' at %d Frequencies',Nfreq);
end
fprintf('\n');

tstart=tic;

for iBlk=1:BlocksNum
    
    fprintf('Process Chn %d to %d. Freq: ',Blocks(iBlk,1),Blocks(iBlk,2));
    
    curChn=Blocks(iBlk,1):Blocks(iBlk,2);
    curChnNum=size(curChn,2);
    
    %set variables in HDR for single channel only
    HDR.InChanSelect=curChn;
    HDR.Calib=HDRin.Calib(1:curChnNum+1,1:curChnNum);
    
    V=sread(HDR,StopSec-StartSec,StartSec); %read whole channel
    
    %demodulate for each frequency
    for iFreq=1:Nfreq
        %display which freq we are doing
        if iFreq < Nfreq
            fprintf('%d,',iFreq);
        else
            fprintf('%d',iFreq);
        end
        % filter and demodulate channel
        [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,Filt{iFreq});
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

