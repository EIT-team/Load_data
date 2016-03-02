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

%find the start time in seconds - subtract 1 to give more samples to reduce
%filtering artefacts
StartSec=floor(InjectionWindows(1,1)/Fs)-1;
if StartSec < 0
    StartSec =0; %make sure it is >=0
end

%find the stop time in seconds - add 1 to give more samples to reduce
%filtering artefacts
StopSec=ceil(InjectionWindows(end,end)/Fs)+1;
if StopSec > Nsec
    StopSec=StopSec-1; %remove added second if we go too far
end


Start_Sample=StartSec*Fs;
Stop_Sample=StopSec*Fs;


if StopSec > Nsec
    error('Stop sec is too long!');
end


%% preallocate

Vmag=nan(size(InjectionWindows,1),Nchn);
Phase=Vmag;




%% Read and Demod each channel
for iChn=1:Nchn
    fprintf('Doing Chn %d of %d\n',iChn,Nchn);
    HDR.InChanSelect=iChn;
    HDR.Calib=HDRin.Calib(1:2,1);
    
    V=sread(HDR,StopSec-StartSec,StartSec); %read whole channel
    [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B,A); % filter and demodulate channel
    [Vmag(:,iChn),Phase(:,iChn)]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod,InjectionWindows-Start_Sample); %process each injection window, adjusting for new start time
    
    
end
disp('done');

%% Pad to complete number of repeats


%find number of repeats - rounding up
Nrep=ceil(size(InjectionWindows,1)/Nprt);

VmagOut=Vmag';
PhaseOut=Phase';

%pad with nans if needed
VmagOut(:,end+1:Nrep*Nprt)=nan;
PhaseOut(:,end+1:Nrep*Nprt)=nan;




%% reshape into correct format
VmagOut=reshape(VmagOut,Nchn*Nprt,Nrep);
PhaseOut=reshape(PhaseOut,Nchn*Nprt,Nrep);


end

