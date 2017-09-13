function [BVstruc] = ScouseTom_ProcessBV( HDR,TT,ExpSetup,BW )
%[BVstruc] = ScouseTom_ProcessBV( HDR,TT,ExpSetup)
%SCOUSETOM_PROCESSBV
%   Demodulates the voltages for "convential" EIT recordings, i.e.
%   sequential pairwise injection, without stimulation. Normally called
%   through ScouseTom_Load
%
%   First the starting line of the injection protocol is estimated, based
%   on the rms of the voltages in the first injection pair. Then the filter
%   settings for all carrier frequencies are determined. 
%
%   The demodulation is then performed on each frequency in turn. The
%   channels are processed in blocks based on an estimate of the total
%   memory usage. The output is then the mean demodulated boundary voltage,
%   as a magnitude BV and Phase PA.
% 
%   The contact impedance is estimated using the voltage on the injection
%   channels and the current amplitude as per the Expsetup.
%
%   The output is stored in a matfile Fname-BV.mat
%
%   Inputs:
%   HDR - from ScouseTom_getHDR contains metadata about EEG file
%   TT  - from ScouseTom_TrigProcess, infomation from digital trigger
%       channels
%   ExpSetup -  Structure containing the setup for the ScouseTom system,
%       this is stored in the FNAME_log.mat with every dataset.
%   BW[100] - Total bandwidth of filter to use i.e. Bandwidth/2 either side
%       of carrier (optional)
%
%   Outputs:
%   BVstruc  - Strucutre containg all boundary voltages, phase angles,
%       contact impedances, along with ExpSetup and filters used furing
%       demodulation. This is also stored in FNAME-BV.mat

%% Do error checking of inputs HERE

%is there any data, do TT and ExpSetup match? Is the freq order ok?
%Expsetup and freq order should match


%% Defaults

if exist('BW','var') == 0  || isempty(BW)
    BW=100; %bandwidth of bandpass filter in demod
end

%% get some variables from inputs Structures

N_elec=size(HDR.InChanSelect,1);
N_freq=ExpSetup.Info.FreqNum;
N_starts=length(TT.InjectionStarts);

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

% some basic info about the data set - HDR is often too big to save each time
info.eegfname=eegfname;
info.TimeNum=datenum(HDR.T0);
info.TimeVec=HDR.T0;

mfilename=fullfile(eegfpath,[eegfname '-BV.mat']);

%create matfile object in same place as data
Out_matfile=matfile(mfilename,'Writable',true);
%% Actually demodulate the data

%find which line in the protocol the data starts with
[StartInj] = ScouseTom_data_checkfirstinj(HDR,TT.InjectionSwitches(1,:),ExpSetup.Protocol );

%find the corresponding filter settings
[Filt,FilterTrim,Fc]=ScouseTom_FindFilterSettings(HDR,TT.InjectionSwitches(1,:),ExpSetup.Protocol(StartInj,1),BW);

%process the data to get the magnitude and phase
[BV,PA,BVSTD,PASTD] = ScouseTom_ReadandDemodChn( HDR,Filt,FilterTrim,TT.InjectionSwitches(1,:),ExpSetup.Protocol,StartInj);

%estimate the contact impedances on the injection electrodes
[Z,Zstd] = ScouseTom_data_estZ( BV,Elec_inj,ZSF);

%% Save the info to the matfile

Out_matfile.TT=TT;% save trigger info
Out_matfile.ExpSetup=ExpSetup; %save the system settings

%save info and the protocol indexes
Out_matfile.keep_idx=keep_idx; % this is the idx of the *measurement* channels which are are interested in
Out_matfile.rem_idx=rem_idx; % this is the idx of *injection* channels, which we are often neglected
Out_matfile.prt_full=prt_full; % the procotol in full [CS+ CS- V+ V-] form, for forward model

%save filter stuff
info.Filt=Filt;
info.trim_demod=FilterTrim;
info.Fc=Fc;
info.ZSF=ZSF;

%save all info
Out_matfile.info=info;

%% save data to matfile

Out_matfile.BV=BV;
Out_matfile.STD=BVSTD;
Out_matfile.Z=Z;
Out_matfile.Zstd=Zstd;
Out_matfile.PhaseAngle=PA;
Out_matfile.PhaseAngleSTD=PASTD;
%% All processing done!

disp('All processing finished! At last!');
teatime=toc(tstart);
fprintf('That took : %.1f seconds \r',teatime);

%output complete structure
BVstruc=load(mfilename);
end

