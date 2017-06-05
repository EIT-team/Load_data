fname='Recording_A.eeg';

% THESE CHANNELS ARE THE ORDER AS ON THE EEG SYSTEM 
Good_chn = [1 3 6 7 10 11 12 13 16 18 19 20 21 23 24 26 28 29 30 31 32];
%Good_chn = [1 3 4 6 7 10 11 12 13 15 16 18 19 20 21 23 24 26 28 29 30 31 32]; %for recording B
%Injection channels, IN ORDER OF FREQUENCIES
Inj_chn = [8,5,9];
Freqs_target=[6000,7000,9000]; % what the current sources were programmed at in arduino DDS code

%%
HDR=sopen(fname,'r',[],['OVERFLOWDETECTION:OFF']);
Fs = HDR.SampleRate;

NumChn=length(Good_chn);
NumFreq = length(Inj_chn);

decimation_factor_eeg =10;
Fs_eeg=Fs/decimation_factor_eeg;
EEG_data=zeros(HDR.SPR/decimation_factor_eeg,NumChn);

decimation_factor_eit =10;
Fs_eit=Fs/decimation_factor_eit;
EIT_data_V=zeros(HDR.SPR/decimation_factor_eit,NumChn*NumFreq);
EIT_data_P=zeros(HDR.SPR/decimation_factor_eit,NumChn*NumFreq);

%% correct for actichamp shifting

for iChn = 1:NumChn
    
    Good_chn_corrected(iChn)=find(strcmp(strtrim(HDR.Label),num2str(Good_chn(iChn))));
    
end

for iChn = 1:length(Inj_chn)
    
    Inj_chn_corrected(iChn)=find(strcmp(strtrim(HDR.Label),num2str(Inj_chn(iChn))));
    
end




%% Filtering Settings
%eeg ones
F_low = 70;
F_high = 0.5;

% eit ones
BW=500; % TOTAL Width
EIT_start_time= 30; % seconds when EIT starts

% find carrier frequencies and then find filter settings for this

for iInj = 1:length(Inj_chn)
    
    HDR=sopen(fname,'r',[Inj_chn_corrected(iInj)],['OVERFLOWDETECTION:OFF']);
    fprintf('Loading data for chn %d...',iInj);
    V=sread(HDR,2,EIT_start_time+5);
    
    %     [ Fc ] = ScouseTom_data_GetCarrier( V,Fs );
    
    [cur_trim_demod,cur_Filt,cur_Fc]=ScouseTom_data_GetFilterTrim(V,Fs,BW,0.03*length(V));
    
    %make it consistent with multifreq bits, which are all cells
    Filt{iInj}=cur_Filt;
    TrimDemod{iInj}=cur_trim_demod;
    Fc{iInj}=cur_Fc;
    
end

%%
tic
for iChn = 1: NumChn
    
    HDR=sopen('Recording_A.eeg','r',[Good_chn_corrected(iChn)],['OVERFLOWDETECTION:OFF']);
    fprintf('Loading data for chn %d...',iChn);
    V=sread(HDR,inf,0);
    
    fprintf('done\n');
    fprintf('Doing EEG processing...');
    %50 Hz notch
    [b,a] = iirnotch(50/(Fs/2),(50/(Fs/2))/45);
    Data = filtfilt(b,a,V);
    % low-pass at F_low
    F_low = 70;
    [b,a] = butter(3,F_low/(Fs/2),'low');
    Data = filtfilt(b,a,Data);
    % High-pass F_high
    F_high = 0.5;
    [b,a] = butter(3,F_high/(Fs/2),'high');
    Data = filtfilt(b,a,Data);
    
    EEG_data(:,iChn)=decimate(Data,decimation_factor_eeg);
    
    fprintf('done\n');
    
    %% EIT FILTERING
    
    fprintf('Doing EIT processing...');
    
    
    for iFreq = 1:NumFreq
        fprintf('%d,',iFreq);
        
        % do the actual demodulation
        [Vdemod,Pdemod]=ScouseTom_data_DemodHilbert(V,Filt{iFreq},TrimDemod{iFreq});
        
        % decimate the voltage and phase signals
        
        Vsig = decimate(Vdemod,decimation_factor_eit);
        Psig = decimate(Pdemod,decimation_factor_eit);
        
        % put this channel at this freq in the correct order
        vidx=(iFreq-1)*NumChn + iChn;
        
        EIT_data_V(:,vidx)=Vsig;
        EIT_data_P(:,vidx)=Psig;
        
    end
    
    fprintf('done\n');
    
    %%
    
    
    
    
    
    toc
end
toc

%% Correct for ActiChamp Gain

% AC gain
AC=load('Freq_AC.mat');
ff=20:20000;
BVdiffFine=spline(AC.F,AC.BVdiff_per_mean,ff); %interpolate data to more freqs
AC_correction=((100-BVdiffFine)/100);


for iFreq = 1:NumFreq
    G1=AC_correction(round(Fc{iFreq}));
    
    vidx=(iFreq-1)*NumChn + 1:(iFreq)*NumChn;
    
    EIT_data_V(:,vidx)=EIT_data_V(:,vidx)*G1;
    
end

save([fname '_processed.mat'],'EIT_data_V','EIT_data_P','EEG_data','Fc','Good_chn','Inj_chn','Freqs_target','TrimDemod','Filt','Fs','Fs_eeg','Fs_eit');

%%
maxtrimsamples = max(cellfun(@(x) max(x),TrimDemod));

BV_mag=mean(Vfull(maxtrimsamples:end-maxtrimsamples,:),1);
BV_phase=mean(Pfull(maxtrimsamples:end-maxtrimsamples,:),1);



T_1 = [1:HDR.SPR/decimation_factor_eeg]/Fs_eeg; %in sec
T = [1:length(Data)]/Fs; %in sec
%%

t = length(Data);
n_t = t/(Fs/100);
for j = 1:10:n_t
    for i = 1:31
        plot(Data(j*Fs:(9+j)*Fs,i) + (i-1)*100);
        hold on;
    end
    hold off;
    waitforbuttonpress;
end