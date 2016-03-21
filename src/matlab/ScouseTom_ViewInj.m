function [  ] = ScouseTom_ViewInj( HDR,ChntoView,Repeat,Startinj,TT,ExpSetup )
%SCOUSETOM_VIEWINJ Summary of this function goes here
%   Detailed explanation goes here


%% check if triggers given

%if TT struct not given then obtain it from HDR
if exist('TT','var') == 0
    Trigger= ScouseTom_TrigReadChn(HDR);
    TT= ScouseTom_TrigProcess(Trigger, HDR);
end

%% check expsetup and injection start

if exist('ExpSetup','var') == 0
   [pathstr,namestr,extstr] = fileparts(HDR.FileName);
mfilename=fullfile(pathstr,[namestr '_log.mat']); 
    %check if _log file is with eeg file
    if exist(mfilename,'file')
        load(mfilename);
    else
        error('Cannot find ExpSetup');
    end
    
end

if exist('StartInj','var') ==0
    [StartInj] = ScouseTom_data_checkfirstinj(HDR,TT.InjectionSwitches(1,:),ExpSetup.Protocol );
end

%% Extract useful stuff
Nelec=size(HDR.InChanSelect,1);


%% Find relevant injection we wish to view

[prt_full,keep_idx,rem_idx,Elec_inj]=ScouseTom_data_findprt(ExpSetup.Protocol,Nelec);

InjtoView=floor(ChntoView/Nelec);
CurrentChannel=ChntoView-(InjtoView*Nelec);




end

