function [ VmagOut,PhaseOut ] = ScouseTom_ReadandDemodChn( HDRin,B,A,Trim_demod,InjectionWindows,StartSec,StopSec,Nprt )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%% Fiddling with inputs
HDR=HDRin;


if strcmp(HDR.TYPE,'BDF');
    %biosemi has status channel as a separate chn in file
    Nchn=HDR.NS -1;
else
    Nchn=HDR.NS;
end

%% Check stuff
Nsec=HDR.NRec;

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
    
    V=sread(HDR,StopSec-StartSec,StartSec);
    [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B,A);
    [Vmag(:,iChn),Phase(:,iChn)]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod,InjectionWindows);
   
   %for loop through each trigger couple, and mean and trim each block
    
    
end
disp('done');

%% other step here?

%need to preallocate and nan pad
%find number of repeats


Nrep=10;
VmagOut=reshape(Vmag',Nchn*Nprt,Nrep);
PhaseOut=reshape(Phase',Nchn*Nprt,Nrep);


end

