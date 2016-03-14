function [ VmagOut,PhaseOut,VmagOutSTD,PhaseOutSTD,Vmag ] = ScouseTom_ReadandDemodChn( HDRin,B,A,Trim_demod,InjectionWindows,Protocol,StartInj )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Fiddling with inputs
HDR=HDRin;

if strcmp(HDR.TYPE,'BDF');
    %biosemi has status channel as a separate chn in file
    %THIS WILL BREAK IF AUX SENSORS SELECTED
    Nchn=HDR.NS -1;
else
    Nchn=HDR.NS;
end

Fs=HDR.SampleRate;

%% Check stuff
Nsec=HDR.NRec;
Nfreq=size(InjectionWindows,2);
Nprt=size(Protocol,1);

if exist('StartInj','var') ==0
    StartInj=ones(Nfreq,1);
    fprintf(2,'WARNING. NO STARTINJ SPECIFIED. ASSUMING 1. \n');
end

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
Stop_Sample=StopSec*Fs;


if StopSec > Nsec
    error('Stop sec is too long!');
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


%% Read and Demod each channel

fprintf('Processing Channels\n');

tstart=tic;

for iChn=1:Nchn
    fprintf('Process Chn: %d of %d. Freq: ',iChn,Nchn);
    
    %set variables in HDR for single channel only
    HDR.InChanSelect=iChn;
    HDR.Calib=HDRin.Calib(1:2,1);
    
    V=sread(HDR,StopSec-StartSec,StartSec); %read whole channel
    
    %demodulate for each frequency
    for iFreq=1:Nfreq
        if iFreq < Nfreq
            fprintf('%d,',iFreq);
        else
            fprintf('%d',iFreq);
        end
        [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B{iFreq},A{iFreq}); % filter and demodulate channel
        [Vmag{iFreq}(:,iChn),PhaseRaw{iFreq}(:,iChn),VmagSTD{iFreq}(:,iChn),PhaseRawSTD{iFreq}(:,iChn)]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod{iFreq},InjectionWindows{iFreq}-Start_Sample); %process each injection window, adjusting for new start time
        
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
    
end






end

