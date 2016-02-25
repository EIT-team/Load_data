function [ output_args ] = ScouseTom_ReadandDemodChn( HDRin,B,A,Trim_demod,InjectionWindows,StartSec,StopSec )
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
Nsec=HDR.Nrec;

if StopSec > NSec
    error('Stop sec is too long!');
end


%% Read and Demod each channel
for iChn=1:Nchn
    
    HDR.NS=iChn;
    
    V=sread(HDR,StartSec,StopSec);
    [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B,A);
    [Vmag,Phase]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod,InjectionWindows);
   
    
   %triggers should be in from [start1, stop1; start2,stop2;...]
   
   %for loop through each trigger couple, and mean and trim each block
    
    
end


%% reshape into each one

%

end

