function [ Filt,TrimDemod,Fc,BaselineCorrection ] = ScouseTom_FindFilterSettings( HDRin,curInjectionSwitches,InjElec,BandWidth )
%SCOUSETOM_FINDFILTERSETTINGS Summary of this function goes here
%   Detailed explanation goes here


%% check if there are actually any injections etc.


if exist('BandWidth','var') ==0
    BandWidth =100;
end


%% Find relevant start inj


Nfreq=size(curInjectionSwitches,2);
Fs=HDRin.SampleRate;

firstinj=zeros(Nfreq,1);
lastinj=zeros(Nfreq,1);
%find earliest switch
for iFreq=1:Nfreq
    firstinj(iFreq)=min(curInjectionSwitches{iFreq}(1,:));
    lastinj(iFreq)=max(curInjectionSwitches{iFreq}(1,:));
end

[FirstSample, FirstFreq]=min(firstinj);
[LastSample, LastFreq]=max(lastinj);

StartSec=floor(FirstSample/Fs);
StopSec=ceil(LastSample/Fs);

StartSample=StartSec*Fs;
StopSample=StopSec*Fs;

%% Load a chunk of data

%preserve original HDR
HDR=HDRin;

%% Find Filter settings for each Freq using relevant data window

for iFreq=1:Nfreq
    
    % load only the injection electrode given
    HDR.InChanSelect=InjElec(iFreq);
    HDR.Calib=HDRin.Calib(1:2,1);
    
    %only load data for first injection
    V=sread(HDR,StopSec-StartSec,StartSec);
    
    %take either the first injection or the first second
    
    tmp=curInjectionSwitches{iFreq}(1,2)-curInjectionSwitches{iFreq}(1,1);
    %     if tmp > Fs
    %         tmp=Fs; %if the first switch is longer than a second, only take a second
    %     end
    
    tmpstart=curInjectionSwitches{iFreq}(1,1)-StartSample;
    tmpidx=tmpstart:tmpstart+tmp;
    
    
    %find carrier frequency and get filter coefficients as well as
    %the amount of data to remove each segment
    [cur_trim_demod,cur_Filt,cur_Fc]=ScouseTom_data_GetFilterTrim(V(tmpidx),Fs,BandWidth );
    
    %make it consistent with multifreq bits, which are all cells
    Filt{iFreq}=cur_Filt;
    TrimDemod{iFreq}=cur_trim_demod;
    Fc{iFreq}=cur_Fc;
end



end

