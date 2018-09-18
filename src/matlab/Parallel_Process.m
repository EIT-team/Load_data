function [EIT,EEG] = Parallel_Process( fname,BWeit)
%PROCESSPARALLEL Summary of this function goes here
%   Detailed explanation goes here

%%
% InjsSim=sort(InjsSim,2);
% Ninj=size(InjsSim,1);

HDR=ScouseTom_getHDR(fname);
Fs=HDR.SampleRate;

if exist('BWeit','var') == 0  || isempty(BWeit)
    BWeit=100; %bandwidth of bandpass filter in demod
end


%% Check number of channels
% get the chan labels from actichamp software
Chn_labels= str2double(HDR.Label);

if any(isnan(Chn_labels))
    Chn_labels=1:size(HDR.Label);
end

% number of channels recorded
Chn_total=max(size(Chn_labels));

% maximum channel number
Chn_max = max(Chn_labels);

gap_in_chn = find((diff(Chn_labels) > 1));

if isempty(gap_in_chn)
    ref_chn = Chn_max +1;
else
    ref_chn = Chn_labels(gap_in_chn) + 1;
end

disp('------------------------------------');
fprintf('Recorded %d channels %d-%d, with reference %d\n',Chn_total,min(Chn_labels),Chn_max,ref_chn);
%% Load the data
StartSec=0;
Secondstoload=inf;

disp('Loading data');
V=sread(HDR,Secondstoload,StartSec);

Vmax= max(V);
Nsamples=size(V,1);


%% Find injection channels

% dont use huge amounts of data for estimation
MaxEstLength=3;
SecondsEst=min([MaxEstLength floor(HDR.SPR/HDR.SampleRate)]);


[Injs, Freqs] = Parallel_FindInjections(V(1:Fs*SecondsEst,:),Fs,Chn_labels);
N_freqs=length(Freqs);


%% Process EIT data
Vfull=nan(size(V,1),N_freqs*Chn_total);
Pfull=Vfull;

BV_inj=nan(Chn_total,N_freqs);
STD_inj=BV_inj;

for iFreq = 1:N_freqs
    fprintf('Processing freq %d\n',iFreq);
    Fc_cur=Freqs(iFreq)+BWeit*[-1 +1];
    
    [cur_Filt,cur_TrimDemod] =ScouseTom_getbpf(20,Fc_cur,Fs);
    
    
    %make it consistent with multifreq bits, which are all cells
    Filt{iFreq}=cur_Filt;
    TrimDemod{iFreq}=cur_TrimDemod;
    Fc{iFreq}=Freqs(iFreq);
    
    
    [Vdemod,Pdemod]=ScouseTom_data_DemodHilbert(V,cur_Filt);
    vidx=(iFreq-1)*Chn_total + 1:(iFreq)*Chn_total;
    
    Vfull(:,vidx)=Vdemod;
    Pfull(:,vidx)=Pdemod;
    
    BV_inj(:,iFreq)=mean(Vdemod(TrimDemod{iFreq}:end-TrimDemod{iFreq},:),1);
    STD_inj(:,iFreq)=std(Vdemod(TrimDemod{iFreq}:end-TrimDemod{iFreq},:),1);
    
end

max_trimsamples = max(cellfun(@(x) max(x),TrimDemod));

Vfull([1:max_trimsamples end-max_trimsamples:end],:)=nan;
Pfull([1:max_trimsamples end-max_trimsamples:end],:)=nan;
%% Correct for ActiChamp Gain

% AC gain
AC=load('Freq_AC.mat');
ff=20:20000;
BVdiffFine=spline(AC.F,AC.BVdiff_per_mean,ff); %interpolate data to more freqs
AC_correction=((100-BVdiffFine)/100);

ACgain=zeros(N_freqs,1);

for iFreq = 1:N_freqs
    G1=AC_correction(round(Fc{iFreq}));
    
    vidx=(iFreq-1)*Chn_total + 1:(iFreq)*Chn_total;
    
    Vfull(:,vidx)=Vfull(:,vidx)*G1;
    
    ACgain(iFreq)=G1;
    
end


disp('Processing Done');


t=(0:length(V)-1)/Fs;

%%

EIT.t=t;
EIT.Vfull=Vfull;
EIT.Pfull=Pfull;
EIT.Injs=Injs;
EIT.Freqs=Freqs;
EIT.info.Fc=Fc;
EIT.info.TrimDemod=TrimDemod;
EIT.info.Filt=Filt;
EIT.info.ACgain=ACgain;



end

