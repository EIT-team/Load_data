function [ output_args ] = ScouseTom_FindFilterSettings( HDRin,TT,InjElec )
%SCOUSETOM_FINDFILTERSETTINGS Summary of this function goes here
%   Detailed explanation goes here


%% check if there are actually any injections etc.


%% only for first injection at the moment


%% Load first second of data

curInjectionSwitches=TT.InjectionSwitches(1);


Nfreq=size(curInjectionSwitches{1,:},2);
Fs=HDR.SampleRate;

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

HDR=HRin;

HDR.InChanSelect=iChn;
HDR.Calib=HDRin.Calib(1:2,1);

%only load data for first injection
V=sread(HDR,StopSec-StartSec,StartSec);

%% find which bit of the dataset to use

%take either the first injection or the first second

tmp=InjectionSwitchesCell{FirstFreq}(1,2)-InjectionSwitchesCell{FirstFreq}(1,1);
if tmp > Fs
    tmp=Fs; %if the first switch is longer than a second, only take a second
end


tmpstart=InjectionSwitchesCell{FirstFreq}(1,1)-StartSample;
tmpidx=tmpstart:tmpstart+tmp;
















tmp=curInjSwitch(idx_f+1)-curInjSwitch(idx_f);
fwind=curInjSwitch(idx_f)-datawindow(1):curInjSwitch(idx_f)-datawindow(1)+tmp;

%find carrier frequency and get filter coefficients as well as
%the amount of data to remove each segment


[trim_demod,B,A,Fc]=ScouseTom_data_GetFilterTrim(V(fwind,Prot(nextprt,1)),Fs,BW,0 );

%make it consistent with multifreq bits, whic are all cells
A={A};
B={B};
trim_demod={trim_demod};




end

