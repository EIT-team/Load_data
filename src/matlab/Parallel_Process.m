function [varargout] = Parallel_Process( fname,BWeit,decimation_factor)
%PROCESSPARALLEL [EIT] = Parallel_Process( fname,BWeit)
%   Processes parallel FDM-EIT data - ASSUMES WHOLE FILE CAN FIT IN RAM
% [EIT] = Parallel_Process(...) % extracts EIT data only
% [EIT,EEG] = Parallel_Process(...) % extracts both EIT and EEG data
%%

% do eeg only if it is asked
if nargout ==2
    DoEEG =1;
else
    DoEEG=0;
end

if exist('BWeit','var') == 0  || isempty(BWeit)
    BWeit=100; %bandwidth of bandpass filter in demod
end

if exist('decimation_factor','var') == 0  || isempty(decimation_factor)
    decimation_factor=1; % default to not doing it, despite the HUGE files this could result in
end

if decimation_factor ==1
    DoDecimation =0;
else
    DoDecimation =1;
end


%% Get Header and triggers

HDR=ScouseTom_getHDR(fname);

if strcmp(HDR.TYPE,'NULL')
    error('File not found');
end

Fs=HDR.SampleRate;

%find triggers
[ StatusChn,TrigPos ] = ScouseTom_geteegtrig(HDR);
% TrigPos = TrigPos/Fs ;


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
t=(0:length(V)-1)/Fs;

Vmax= max(V);
Nsamples=size(V,1);
%% Find injection channels

% dont use huge amounts of data for estimation
MaxEstLength=3;
SecondsEst=min([MaxEstLength floor(HDR.SPR/HDR.SampleRate)]);

[Injs, Freqs] = Parallel_FindInjections(V(1:Fs*SecondsEst,:),Fs,Chn_labels);
N_freqs=length(Freqs);
%% Find Decimation factors

if DoDecimation
    %decimation of 13 or more should be done by separate passes.
    if decimation_factor > 13
        
        facs= factor(decimation_factor);
        
        % if factor is bigger than 13 then this doesnt work, so better not to
        % do anything
        if any(facs > 13)
            error('Too large prime factor! Adjust so all are smaller than 13');
        end
        
        % split the decimation factor into steps less than 12
        next_div=inf;
        
        decimation_factor_vec=[1];
        vec_cnt=1;
        while next_div > 1
            
            %find the next biggest divisor thats less than 12
            cur_divisors=divisors(decimation_factor/(prod([decimation_factor_vec])));
            next_div=cur_divisors(find(cur_divisors < 12,1,'last'));
            
            if next_div > 1
                decimation_factor_vec(vec_cnt)=next_div;
                vec_cnt=vec_cnt+1;
            end
        end
        
        
    else
        decimation_factor_vec = decimation_factor;
    end
end
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

%% Process EEG data

if DoEEG
    disp('Processing EEG')
    
    [EEGlpFilt, TrimDemodEEGlpf] = ScouseTom_getlpf(6,400,Fs);
    [EEGhpFilt, TrimDemodEEGhpf] = ScouseTom_gethpf(1,2,Fs);
    
    EEG_data=filtfilt(EEGlpFilt,V);
    EEG_data=filtfilt(EEGhpFilt,EEG_data);
    
end

%% Decimate

if DoDecimation
    fprintf('Decimating...');
    
    Fs_dec=Fs/decimation_factor;
    Nsamples_dec=Nsamples/decimation_factor;
    
    fprintf('EIT...');
    
    V_dec=nan(Nsamples_dec,size(Vfull,2));
    P_dec=nan(Nsamples_dec,size(Vfull,2));
    
    for iChn = 1:size(Vfull,2)
        Vtmp=Vfull(:,iChn);
        Ptmp=Pfull(:,iChn);
        for iDec = 1:length(decimation_factor_vec)
            Vtmp=decimate(Vtmp,decimation_factor_vec(iDec));
            Ptmp=decimate(Ptmp,decimation_factor_vec(iDec));
        end
        V_dec(:,iChn)=Vtmp;
        P_dec(:,iChn)=Ptmp;
    end
    
    if DoEEG
        fprintf('EEG...');
        EEG_data_dec=nan(Nsamples_dec,Chn_max);
        for iChn = 1:Chn_max
            Vtmp=EEG_data(:,iChn);
            for iDec = 1:length(decimation_factor_vec)
                Vtmp=decimate(Vtmp,decimation_factor_vec(iDec));
            end
            EEG_data_dec(:,iChn)=Vtmp;
        end
    end
    
    
    % Do triggers as well
    
    
    
    
    
    fprintf('done\n');
end
%% Trim data
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
%% Find protocol info
[prt_full,keep_idx,rem_idx]=ScouseTom_data_findprt(Injs,Chn_total); %from ScouseTom Repo




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
EIT.info.TrimMax=max_trimsamples;
EIT.protocol.prt=prt_full;
EIT.protocol.keep_idx=keep_idx;
EIT.protocol.rem_idx=rem_idx;
EIT.Trig.TrigPos=TrigPos;


EEG.t=t;
EEG.Data=EEG_data;
EEG.info.EEGlpf=EEGlpFilt;
EEG.info.EEGhpf=EEGhpFilt;
EEG.info.TrimDemodEEGlpf = TrimDemodEEGlpf;
EEG.info.TrimDemodEEGhpf = TrimDemodEEGhpf;
EEG.Trig.TrigPos=TrigPos;

varargout{1}=EIT;
if DoEEG
    varargout{2}=EEG;
end

end

