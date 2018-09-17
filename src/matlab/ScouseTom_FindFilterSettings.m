function [ Filt,TrimDemod,Fc ] = ScouseTom_FindFilterSettings( HDRin,curInjectionSwitches,InjElec,BandWidth )
% ScouseTom_FindFilterSettings( HDRin,curInjectionSwitches,InjElec,BandWidth )
%SCOUSETOM_FINDFILTERSETTINGS Finds the optimal filters for all injected
%frequencies in a given file. This is based on the carrier frequency and
%the length of injection. If possible, an IIR filter with smaller number of
%coefficients will be chosen. However if the injections are too short, and
%the filter has not sufficiently decayed, an FIR filter will be used
%instead.
%
% Most common usage is like this:
% [Filt,FilterTrim,Fc]=ScouseTom_FindFilterSettings(HDR,TT.InjectionSwitches(1,:),ExpSetup.Protocol(StartInj,1),BW);
%
% This function is basically a wrapper for ScouseTom_data_GetFilterTrim
% which handles all of the filter settings.
%
%
% Inputs:
% HDRin - HDR structure from ScouseTom_getHDR
%
% curInjectionSwitches - the timepoints (in samples) which correspond the
%   swtiching of injection pairs in the data. This is found in
%   TT.InjectionSwitches, the output of ScouseTom_TrigProcess
%
% InjElec - vector with the electrode numbers corresponding to the first
%   injection pair i.e. [1 16]. We need this as these are the largest
%   channels and thus best suited for finding the carrier freq
%
% BandWidth - Total bandwidth of filter to use i.e. Bandwidth/2 either side
% of carrier
%
% Outputs:
% Filt - filter object (used like filtfilt(Filt,V)). This is a
%   cell if multiple frequencies in file.
% TrimDemod - The number of samples to remove before and after a switch of injection pair before averaging.
%   To remove the errors from filtering artefacts. A cell if mutiple
%   frequencies used
% Fc - Carrier frequencies, cell if multiple freqs used


%% set default

% This is TOTAL bandwidth, i.e. BandWidth/2 either side of the carrier
if exist('BandWidth','var') ==0
    BandWidth =100;
end


%% Find relevant start injection

%number of frequencies
Nfreq=size(curInjectionSwitches,2);
Fs=HDRin.SampleRate;

firstinj=zeros(Nfreq,1);
lastinj=zeros(Nfreq,1);
%find earliest switch
for iFreq=1:Nfreq
    firstinj(iFreq)=min(curInjectionSwitches{iFreq}(1,:));
    lastinj(iFreq)=max(curInjectionSwitches{iFreq}(1,:));
end

% find the chunk of data we want to load
[FirstSample, FirstFreq]=min(firstinj);
[LastSample, LastFreq]=max(lastinj);

%round to nearest whole second (this is what sread needs)
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
    
    %only load data for first injection rounded to nearest seconds
    V=sread(HDR,StopSec-StartSec,StartSec);
    
    %take the first injection switch indicies
    tmp=curInjectionSwitches{iFreq}(1,2)-curInjectionSwitches{iFreq}(1,1);
    
    % injection switches are with respect to whole file, so adjust for the
    % start of the data we have collected
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

