function [BVstruc] = ScouseTom_ProcessBV( HDR,TT,ExpSetup,varargin )
%SCOUSETOM_ Summary of this function goes here
%   Detailed explanation goes here

% ONLY DOES ONE INJECTION AT THE MOMENT

%% Do error checking of inputs HERE

%is there any data, do TT and ExpSetup match? Is the freq order ok?
%Expsetup and freq order should match

%allow for passing filter info into function


%% Defaults
BW=50; %bandwidth of bandpass filter in demod


%% See what type of system we are using


%the maximum voltage is different for each system, this *should* be in the
%HDR structure somewhere, but I dont know where it is

switch HDR.TYPE
    case 'BDF' % biosemi file
        MaxV=0.5e6; %500mV range on BioSemi
    case 'EEG'
end


%% get some variables from inputs Structures

Prot=ExpSetup.Protocol;
N_prt=size(Prot,1);
N_elec=ExpSetup.Elec_num;
N_freq=ExpSetup.Info.FreqNum;
N_starts=length(TT.InjectionStarts);

Fs=HDR.SampleRate;
eegfname=HDR.FILE.Name;
eegfpath=HDR.FILE.Path;


fprintf('Processing data in %s\n',eegfname);
tstart=tic;



%% calculate the keep and rem idx
% need this beore loading data to know which channels to estimate contact
% impedance on

[prt_full,keep_idx,rem_idx,Elec_inj]=ScouseTom_data_findprt(ExpSetup.Protocol,N_elec);

%scale factor - impedance conversion
ZSF=1./(ExpSetup.Amp); %keep this in uA as voltages are in uV

%% create matfile object for saving data

%matfile object can handle big files, and stops everything being stored in
%memory

% some basic info about the data set - HDR is too big to save each time
info.eegfname=eegfname;
info.TimeNum=datenum(HDR.T0);
info.TimeVec=HDR.T0;

%create matfile object in same place as data
bigmat=matfile(fullfile(eegfpath,[eegfname '-BV.mat']),'Writable',true);


%% Loop through each injection start in file

%find which line in the protocol the data starts with
[StartInj] = ScouseTom_data_checkfirstinj(HDR,TT.InjectionSwitches(1,:),ExpSetup.Protocol );

%find the corresponding filter settings
[B,A,FilterTrim,Fc]=ScouseTom_FindFilterSettings(HDR,TT.InjectionSwitches(1,:),ExpSetup.Protocol(StartInj,1));

%process the data to get the magnitude and phase
[BV,PA,BVSTD,PASTD] = ScouseTom_ReadandDemodChn( HDR,B,A,FilterTrim,TT.InjectionSwitches(1,:),ExpSetup.Protocol);

%estimate the contact impedances on the injection electrodes
[Z,Zstd] = ScouseTom_data_estZ( BV,Elec_inj,ZSF);



%% Save the info to the matfile

% save trigger info
bigmat.TT=TT;
%save the system settings
bigmat.ExpSetup=ExpSetup;
%save info and the protocol indexes
bigmat.keep_idx=keep_idx;
bigmat.rem_idx=rem_idx;
bigmat.prt_full=prt_full;

%save filter stuff
info.B=B;
info.A=A;
info.trim_demod=FilterTrim;
info.Fc=Fc;
info.ZSF=ZSF;

%save all info
bigmat.info=info;

%% save data to matfile

bigmat.BV=BV;
bigmat.STD=BVSTD;
bigmat.Z=Z;
bigmat.Zstd=Zstd;
bigmat.PhaseAngle=PA;
bigmat.PhaseAngleSTD=PASTD;
%% All processing done!

disp('All processing finished! At fucking last!');
teatime=toc(tstart);
fprintf('That took : %.1f seconds \r',teatime);

%output complete structure
BVstruc=bigmat.BV;
end

