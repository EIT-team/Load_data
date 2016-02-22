function [ output_args ] = ScouseTom_ReadandDemodChn( HDRin,B,A,StartSec,StopSec )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


HDR=HDRin;


if strcmp(HDR.TYPE,'BDF');
    %biosemi has status channel as a separate chn in file
Nchn=HDR.NS -1;
else
    Nchn=HDR.NS;
end


Nsec=HDR.Nrec;

if StopSec > NSec
    error('Stop sec is too long!');
end

for iChn=1:Nchn
    
    HDR.NS=iChn;
    
    V=sread(HDR,StartSec,StopSec);
    







end

