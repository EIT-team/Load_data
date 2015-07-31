function [ output_args ] = ScouseTom_ProcessBV( HDR,TT,ExpSetup )
%SCOUSETOM_ Summary of this function goes here
%   Detailed explanation goes here


%% Do error checking of inputs HERE

%is there any data, do TT and ExpSetup match? 



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

Fs=HDR.SampleRate;

eegfname=HDR.FILE.Name;
eegfpath=HDR.FILE.Path;


%% calculate the keep and rem idx 
% need this beore loading data to know which channels to estimate contact
% impedance on

injs=ExpSetup.Protocol; %injection pairs
chn=N_elec; %number of channels
vp=(1:chn)';% positive voltage channel 1  to Number of electrodes
vm=ones(size(vp))*(chn+1); %always against a ground electrode

%make the entire protocol each line is INJ+ INJ- MEAS+ MEAS-
prt_mat=[];
for iii=1:size(injs,1)
    temp=[repmat(injs(iii,:),chn,1) vp vm];
    prt_mat=[prt_mat ; temp];
end

%find remove index, any protocol lines including the injetion channels are
%"bad"
prt_full=prt_mat;
prt=prt_full;
rem_idx=[];
for iPrt = 1:size(prt,1)
    if any(ismember(prt_full(iPrt,1:2),prt(iPrt,3:4))) ==1
        rem_idx=[rem_idx,iPrt];
    end
end
%keep index is anything that we *dont* remove
keep_idx=setdiff(1:length(prt_full),rem_idx);

%% get injection channels for use in contact impedance calculations

%loop through protocol - find which lines in the BV the injection
%electrodes belong to then add then to an array for each electrode. 

%I cant remember why I do this separately to the stuff about injections
%above...

%Electrode Injections
Elec_inj=nan(N_elec,N_prt);

for iPrt = 1:N_prt
    Prt_cur=Prot(iPrt,:);
    start_idx=((iPrt-1)*N_elec);
    BV_chn=start_idx+Prt_cur;
    Elec_inj(Prt_cur,iPrt)=BV_chn;
end

Elec_inj=sort(Elec_inj,2);

%clear up matrix
Elec_inj(:,all(isnan(Elec_inj),1))=[];

%scale factor - impedance conversion
ZSF=1/((1e6)*ExpSetup.Amp);

%% create matfile object for saving data THIS IS WHERE I STOPPED

info.eegfname=eegfname;
info.TimeNum=datenum(HDR.T0);
info.TimeVec=HDR.T0;

bigmat=matfile(fullfile(eegfpath,[eegfname '-BV.mat']),'Writable',true);

bigmat.SWW=SWW;

bigmat.ExpSetup=ExpSetup;

bigmat.info=info;

bigmat.keep_idx=keep_idx;
bigmat.rem_idx=rem_idx;
bigmat.prt_full=prt_full;





%% Start loading stuff




end

