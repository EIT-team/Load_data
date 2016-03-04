function [ VmagOut,PhaseOut ] = ScouseTom_ReadandDemodChn( HDRin,B,A,Trim_demod,InjectionWindows,Nprt )
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
Phase=Vmag;




%% Read and Demod each channel
for iChn=1:Nchn
    fprintf('Doing Chn %d of %d\n',iChn,Nchn);
    HDR.InChanSelect=iChn;
    HDR.Calib=HDRin.Calib(1:2,1);
    
    V=sread(HDR,StopSec-StartSec,StartSec); %read whole channel
    
    %demodulate for each frequency
    for iFreq=1:Nfreq
        
        
        
        
        
        
        
        
        [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B{iFreq},A{iFreq}); % filter and demodulate channel
        [Vmag{iFreq}(:,iChn),Phase{iFreq}(:,iChn)]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod{iFreq},InjectionWindows{iFreq}-Start_Sample); %process each injection window, adjusting for new start time
        
    end
    
end
disp('done');

%% Pad to complete number of repeats

for iFreq=1:Nfreq
    
    %find number of repeats - rounding up
    Nrep=ceil(size(InjectionWindows{iFreq},1)/Nprt);
    
    VmagOutTmp=Vmag{iFreq}';
    PhaseOutTmp=Phase{iFreq}';
    
    %pad with nans if needed
    VmagOutTmp(:,end+1:Nrep*Nprt)=nan;
    PhaseOutTmp(:,end+1:Nrep*Nprt)=nan;
    
    %% reshape into correct format
    VmagOut{iFreq}=reshape(VmagOutTmp,Nchn*Nprt,Nrep);
    PhaseOut{iFreq}=reshape(PhaseOutTmp,Nchn*Nprt,Nrep);
    
end






end

