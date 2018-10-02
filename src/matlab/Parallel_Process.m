function [varargout] = Parallel_Process( fname,BWeit,decimation_factor)
%PROCESSPARALLEL [EIT] = Parallel_Process( fname,BWeit)
%   Processes parallel FDM-EIT data - ASSUMES WHOLE FILE CAN FIT IN RAM
% [EIT] = Parallel_Process(...) % extracts EIT data only
% [EIT,EEG] = Parallel_Process(...) % extracts both EIT and EEG data
%% Process Inputs

% do eeg only if it is asked
if nargout ==2
    DoEEG =1;
else
    DoEEG=0;
end

if BWeit ==0 % dont do EIT - when EIT not present, things get confused when you end up with 50Hz injections
    DoEIT =0;
else
    DoEIT =1;
end
%set bandwidth
if exist('BWeit','var') == 0  || isempty(BWeit)
    BWeit=100; %bandwidth of bandpass filter in demod
end

% dont do decimation by default
if exist('decimation_factor','var') == 0  || isempty(decimation_factor) || decimation_factor==0
    decimation_factor=1; % default to not doing it, despite the HUGE files this could result in
end

%choose decimation factor
if decimation_factor ==1
    DoDecimation =0;
else
    DoDecimation =1;
    
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
%% Get Header and triggers

if ischar(fname)
    
    HDR=ScouseTom_getHDR(fname);
else if isstruct(fname)
        HDR=fname;
    end
end

if strcmp(HDR.TYPE,'NULL')
    error('File not found');
end

Fs=HDR.SampleRate;

%find triggers
[ StatusChn,TrigPos ] = ScouseTom_geteegtrig(HDR);
% TrigPos = TrigPos/Fs ;

% check for non integer decimation factor here so it breaks early
Fs_dec=Fs/decimation_factor;

if Fs_dec ~= int32(Fs_dec)
    error(sprintf('Non integer decimated sample rate, Fs:%d Fs_dec:%f',Fs,Fs_dec));
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

%% Process EIT data

if DoEIT
    
    %%  Find injection channels
    % dont use huge amounts of data for estimation
    MaxEstLength=3;
    SecondsEst=min([MaxEstLength floor(HDR.SPR/HDR.SampleRate)]);
    
    [Injs, Freqs] = Parallel_FindInjections(V(1:Fs*SecondsEst,:),Fs,Chn_labels);
    N_freqs=length(Freqs);
    %% process each frequency in turn
    %preallocate
    EIT_data_V=nan(size(V,1),N_freqs*Chn_total);
    EIT_data_P=EIT_data_V;
    Filt=cell(1,N_freqs);
    TrimDemod=zeros(1,N_freqs);
    Fc=TrimDemod;
    
    for iFreq = 1:N_freqs
        fprintf('Processing freq %d\n',iFreq);
        Fc_cur=Freqs(iFreq)+BWeit*[-1 +1];
        
        [cur_Filt,cur_TrimDemod] =ScouseTom_getbpf(20,Fc_cur,Fs);
        
        %make it consistent with multifreq bits, which are all cells
        Filt{iFreq}=cur_Filt;
        TrimDemod(iFreq)=cur_TrimDemod;
        Fc(iFreq)=Freqs(iFreq);
        
        [Vdemod,Pdemod]=ScouseTom_data_DemodHilbert(V,cur_Filt);
        vidx=(iFreq-1)*Chn_total + 1:(iFreq)*Chn_total;
        
        EIT_data_V(:,vidx)=Vdemod;
        EIT_data_P(:,vidx)=Pdemod;
    end
    %% Correct for ActiChamp Gain
    
    % AC gain
    AC=load('Freq_AC.mat');
    ff=20:20000;
    BVdiffFine=spline(AC.F,AC.BVdiff_per_mean,ff); %interpolate data to more freqs
    AC_correction=((100-BVdiffFine)/100);
    
    ACgain=zeros(N_freqs,1);
    
    for iFreq = 1:N_freqs
        G1=AC_correction(round(Fc(iFreq)));
        
        vidx=(iFreq-1)*Chn_total + 1:(iFreq)*Chn_total;
        
        EIT_data_V(:,vidx)=EIT_data_V(:,vidx)*G1;
        
        ACgain(iFreq)=G1;
    end
    
    %% Find protocol info
    [prt_full,keep_idx,rem_idx]=ScouseTom_data_findprt(Injs,Chn_total); %from ScouseTom Repo
    
else
    N_freqs=1;
    TrimDemod=0;
end

%% Process EEG data

if DoEEG
    disp('Processing EEG')
    
    [EEGlpFilt, TrimDemodEEGlpf] = ScouseTom_getlpf(6,400,Fs);
    [EEGhpFilt, TrimDemodEEGhpf] = ScouseTom_gethpf(1,2,Fs);
    
    
    TrimDemod=[TrimDemod TrimDemodEEGlpf TrimDemodEEGhpf ];
    
    EEG_data=filtfilt(EEGlpFilt,V);
    EEG_data=filtfilt(EEGhpFilt,EEG_data);
    
end
%% Decimate

if DoDecimation
    fprintf('Decimating...');
    Nsamples_dec=Nsamples/decimation_factor;
    if DoEIT
        
        fprintf('EIT...');
        
        V_dec=nan(Nsamples_dec,size(EIT_data_V,2));
        P_dec=nan(Nsamples_dec,size(EIT_data_V,2));
        
        for iChn = 1:size(EIT_data_V,2)
            Vtmp=EIT_data_V(:,iChn);
            Ptmp=EIT_data_P(:,iChn);
            for iDec = 1:length(decimation_factor_vec)
                Vtmp=decimate(Vtmp,decimation_factor_vec(iDec),100,'fir');
                Ptmp=decimate(Ptmp,decimation_factor_vec(iDec),100,'fir');
            end
            V_dec(:,iChn)=Vtmp;
            P_dec(:,iChn)=Ptmp;
        end
        %replace variables with decimated ones
        EIT_data_V=V_dec;
        EIT_data_P=P_dec;
    end
    
    if DoEEG
        fprintf('EEG...');
        EEG_data_dec=nan(Nsamples_dec,Chn_max);
        for iChn = 1:Chn_max
            Vtmp=EEG_data(:,iChn);
            for iDec = 1:length(decimation_factor_vec)
                Vtmp=decimate(Vtmp,decimation_factor_vec(iDec),100,'fir');
            end
            EEG_data_dec(:,iChn)=Vtmp;
        end
        %replace variables with decimated ones
        EEG_data=EEG_data_dec;
    end
    % Do triggers as well
    TrigPos=round(TrigPos/decimation_factor);
    
    %replace variables with decimated ones
    Fs=Fs_dec;
    TrimDemod=TrimDemod/decimation_factor;
    Nsamples=Nsamples_dec;
    clear V_dec P_dec EEG_data_dec
    
    fprintf('done\n');
end
%% Trim data

max_trimsamples = max(TrimDemod);
max_trimsamples = ceil(max_trimsamples/decimation_factor); % correct for decimation factor
if DoEIT
    EIT_data_V([1:max_trimsamples end-max_trimsamples:end],:)=nan;
    EIT_data_P([1:max_trimsamples end-max_trimsamples:end],:)=nan;
end
if DoEEG
    EEG_data([1:max_trimsamples end-max_trimsamples:end],:)=nan;
end
disp('Processing Done');
%% Save everything to Strucutre

t=(0:Nsamples-1)/Fs;

if DoEIT
    EIT.t=t;
    EIT.Fs=Fs;
    EIT.Data_V=EIT_data_V;
    EIT.Data_P=EIT_data_P;
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
else
    EIT=nan;
end
varargout{1}=EIT;


if DoEEG
    
    EEG.t=t;
    EEG.Fs=Fs;
    EEG.Data=EEG_data;
    EEG.info.EEGlpf=EEGlpFilt;
    EEG.info.EEGhpf=EEGhpFilt;
    EEG.info.TrimDemodEEGlpf = TrimDemodEEGlpf;
    EEG.info.TrimDemodEEGhpf = TrimDemodEEGhpf;
    EEG.Trig.TrigPos=TrigPos;
    
    varargout{2}=EEG;
end

end

