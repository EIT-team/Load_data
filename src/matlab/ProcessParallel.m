function [EEG_dV_signal,EIT_dV_signal, t_signal,EIT_signal_V,EIT_baseline_V,EEG_data,Vmax] = ProcessParallel( fname,BaselineWindow,SignalWindow,TimeStep,BWeit,InjsSim,BV0,F,plotflag)
%PROCESSPARALLEL Summary of this function goes here
%   Detailed explanation goes here

%%
InjsSim=sort(InjsSim,2);
Ninj=size(InjsSim,1);

if exist('plotflag','var') == 0
    plotflag = 1;
end

HDR=ScouseTom_getHDR(fname);

Fs=HDR.SampleRate;
% Nchn=size(HDR.InChanSelect,1);
Nchn = size(cell2mat(regexp(strtrim(HDR.Label),'Ch')),1); %this is different if its recorded by pycorder to activew...

[ StatusChn,TrigPos ] = ScouseTom_geteegtrig(HDR);
TrigPos = TrigPos/Fs ;

Trig_gaps = [ round( diff(TrigPos))];

if isempty(BaselineWindow)
    %if not specified use the first big gap in the triggers
    BaselineWindow(1) = (TrigPos(find (Trig_gaps > 1, 1)));
    BaselineWindow(2) = (TrigPos(find (Trig_gaps > 1, 1)+1));
    
    fprintf('Taking baseline between %.1f and %.1f s\n',BaselineWindow(1),BaselineWindow(2));
    
end

if isempty(SignalWindow)
    %if not specified use the first big gap in the triggers
    SignalWindow(1) = (TrigPos(find (Trig_gaps > 1, 1,'last')));
    SignalWindow(2) = (TrigPos(find (Trig_gaps > 1, 1,'last')+1));
    fprintf('Taking signal between %.1f and %.1f s\n',SignalWindow(1),SignalWindow(2));
end


% StartSec=max([min([floor(BaselineWindow) floor(SignalWindow)])-4 0]);
StartSec=0;
StopSec=max([ceil(BaselineWindow) ceil(SignalWindow)])+4;
StopSec =min([StopSec ( floor(HDR.SPR/HDR.SampleRate)-1)]);

%% EIT Filter Bandwidth and Decimation factors

% sample rate we want out based on time steps
Fs_target=1/TimeStep;

% filter bandwidth this corresponds to during demodulation
BW = 2 * Fs_target;
BW= BWeit;  % use the one given by user instead as we have to have more time steps in order to capture the EEG signal, this would give a too wide BW for the EIT
decimation_factor =Fs/Fs_target;

% find the decimation factors
if decimation_factor > 13
    
    facs= factor(decimation_factor);
    
    if any(facs > 13)
        error('Too large prime factor! Adjust so all are smaller than 13');
    end
    
    % find some way of reducing these
    %%
    next_div=inf;
    
    decimation_factor_vec=[1];
    
    vec_cnt=1;
    
    while next_div > 1
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



%% EEG filters

lpFilt = designfilt('lowpassiir', ...       % Response type
    'PassbandFrequency',500, ...
    'StopbandFrequency',800, ...
    'StopbandAttenuation',120, ...   % Magnitude constraints
    'PassbandRipple',0.1, ...
    'DesignMethod','butter', ...      % Design method
    'MatchExactly','passband', ...   % Design method options
    'SampleRate',Fs)     ;          % Sample rate

hpFilt = designfilt('highpassiir', ...
    'PassbandFrequency',2, ...
    'PassbandRipple',1, ...
    'StopbandFrequency',0.1, ...
    'StopbandAttenuation',30, ...   % Magnitude constraints
    'DesignMethod','butter', ...      % Design method
    'MatchExactly','stopband', ...   % Design method options
    'SampleRate',Fs);



%% Load the data
disp('loading data');
V=sread(HDR,StopSec-StartSec,StartSec);

Vmax= max(V);

if plotflag
    figure
    bar(Vmax)
    xlabel('channel')
    ylabel('uV')
    title('max voltage - USE THIS TO SEE WHAT ELECS YOU WANT TO REMOVE!!!');
    ylim(max(ylim) * [-1 1]);
    drawnow
end

Nsamples=size(V,1);

%% Time chunks
baseline_wind=BaselineWindow-StartSec;
signal_wind=SignalWindow-StartSec;

%%
Fs_dec=Fs/decimation_factor;
EEG_data=zeros(Nsamples/decimation_factor,Nchn);


EIT_data_V=zeros(Nsamples/decimation_factor,Nchn*Ninj);
EIT_data_P=zeros(Nsamples/decimation_factor,Nchn*Ninj);


%% Find injection electrodes and freqs

if exist('F','var') == 0
    
    [InjsExp, Freqs] = Find_Injection_Freqs_And_Elecs(V(1:Fs,:),Fs);
    
    %make sure they are in the same order
    InjsExp=sort(InjsExp,2);
    InjsSim=sort(InjsSim,2);
    
    %put freqs in order to give same protocol as simulation
    [~,Locb]=ismember(InjsSim(:,1),InjsExp(:,1));
    F=Freqs(Locb);
    
end

nFreq=length(F);

%% Filtering Settings

% find carrier frequencies and then find filter settings for this

if ~(BWeit==0)
    
    
    for iInj = 1:Ninj
        
        [cur_trim_demod,cur_Filt,cur_Fc]=ScouseTom_data_GetFilterTrim(V(1:Fs,InjsSim(iInj,1)),Fs,BW,3*Fs);
        
        %make it consistent with multifreq bits, which are all cells
        Filt{iInj}=cur_Filt;
        TrimDemod{iInj}=cur_trim_demod;
        Fc{iInj}=cur_Fc;
        
    end
    
end

%% do EEG Stuff

fprintf('Doing EEG processing...');
% Data_eeg_filt = filtfilt(hpFilt,V);
% Data_eeg_filt = filtfilt(lpFilt,Data_eeg_filt);
% %%
% fprintf('Decimating...');
% for iChn = 1:Nchn
%     Vtmp=Data_eeg_filt(:,iChn);
%     for iDec = 1:length(decimation_factor_vec)
%         Vtmp=decimate(Vtmp,decimation_factor_vec(iDec));
%     end
%     EEG_data(:,iChn)=Vtmp;
% end
fprintf('done\n');


% clear Data_eeg_filt
%% Do EIT stuff
if ~(BWeit==0)
    
    fprintf('Doing EIT processing...');
    
    for iFreq = 1:Ninj
        fprintf('%d,',iFreq);
        
        % do the actual demodulation
        [Vdemod,Pdemod]=ScouseTom_data_DemodHilbert(V,Filt{iFreq},TrimDemod{iFreq});
        
        % decimate the voltage and phase signals
        
        for iChn = 1:Nchn
            Vtmp=Vdemod(:,iChn);
            Ptmp=Pdemod(:,iChn);
            for iDec = 1:length(decimation_factor_vec)
                Vtmp=decimate(Vtmp,decimation_factor_vec(iDec),100,'fir');
                Ptmp=decimate(Ptmp,decimation_factor_vec(iDec),100,'fir');
            end
            
            % put this channel at this freq in the correct order
            vidx=(iFreq-1)*Nchn + iChn;
            EIT_data_V(:,vidx)=Vtmp;
            EIT_data_P(:,vidx)=Ptmp;
        end
        
    end
    
    fprintf('done\n');
    
    %% Correct for ActiChamp Gain
    
    % AC gain
    AC=load('Freq_AC.mat');
    ff=20:20000;
    BVdiffFine=spline(AC.F,AC.BVdiff_per_mean,ff); %interpolate data to more freqs
    AC_correction=((100-BVdiffFine)/100);
    
    
    for iFreq = 1:Ninj
        G1=AC_correction(round(Fc{iFreq}));
        
        vidx=(iFreq-1)*Nchn + 1:(iFreq)*Nchn;
        
        EIT_data_V(:,vidx)=EIT_data_V(:,vidx)*G1;
        
    end
    
    
    
end
% t=(0:length(V)-1)/Fs;
t_1=(0:(length(V))/decimation_factor-1)/Fs_dec;




%% Correct units

EIT_data_V=1e-6*EIT_data_V.*(sign(BV0'));

EEG_data_V=1e-6*EEG_data;



%% plot data here
if plotflag
    figure;
    hold on
    h1=plot(StartSec+t_1,EIT_data_V);
    h2=plot(StartSec+[baseline_wind(1) baseline_wind(1)],ylim,'k--','DisplayName','BaselineWindow');
    h3=plot(StartSec+[baseline_wind(2) baseline_wind(2)],ylim,'k--');
    h4=plot(StartSec+[signal_wind(1) signal_wind(1)],ylim,'r:','DisplayName','SignalWindow');
    h5=plot(StartSec+[signal_wind(2) signal_wind(2)],ylim,'r:');
    hold off
    legend([h2 h4])
    title('Demodulated signals')
    drawnow
end


%% take chunks

EIT_baseline_V=mean(EIT_data_V(t_1 > baseline_wind(1) & t_1 < baseline_wind(2),:));
EIT_baseline_P=mean(EIT_data_P(t_1 > baseline_wind(1) & t_1 < baseline_wind(2),:));
EEG_baseline=mean(EEG_data(t_1 > baseline_wind(1) & t_1 < baseline_wind(2),:));
% P_baseline=mean(Pbin(tbin > baseline_wind(1) & tbin < baseline_wind(2),:));

signal_idx=t_1 > signal_wind(1) & t_1 < signal_wind(2);

EIT_signal_V=EIT_data_V(signal_idx,:);
EIT_signal_P=EIT_data_P(signal_idx,:);


EEG_signal=EEG_data_V(signal_idx,:);




t_signal=(0:size(signal_idx)-1)/Fs;

%% change in voltage

EIT_dV_full=(EIT_data_V) - EIT_baseline_V;
EIT_dV_signal=EIT_dV_full(signal_idx,:);

EEG_dV_full=(EEG_data) - EEG_baseline;
EEG_dV_signal=EEG_dV_full(signal_idx,:);


%% plot final dV
if plotflag
    [maxval, maxchn]=max(max(EIT_dV_signal));
    maxval = max([ 1 maxval]);
    
    figure
    hold on
    h1=plot(StartSec+t_1,EIT_dV_full);
    ylim(maxval*[-1 1])
    h2=plot(StartSec+[baseline_wind(1) baseline_wind(1)],ylim,'k--','DisplayName','BaselineWindow');
    h3=plot(StartSec+[baseline_wind(2) baseline_wind(2)],ylim,'k--');
    h4=plot(StartSec+[signal_wind(1) signal_wind(1)],ylim,'r:','DisplayName','SignalWindow');
    h5=plot(StartSec+[signal_wind(2) signal_wind(2)],ylim,'r:');
    hold off
    legend([h2 h4])
    title('EIT Voltage Change whole data set')
    drawnow
    
    [maxval, maxchn]=max(max(EEG_dV_signal));
    
    figure
    hold on
    h1=plot(StartSec+t_1,EEG_dV_full);
%     ylim(maxval*[-1 1])
    h2=plot(StartSec+[baseline_wind(1) baseline_wind(1)],ylim,'k--','DisplayName','BaselineWindow');
    h3=plot(StartSec+[baseline_wind(2) baseline_wind(2)],ylim,'k--');
    h4=plot(StartSec+[signal_wind(1) signal_wind(1)],ylim,'r:','DisplayName','SignalWindow');
    h5=plot(StartSec+[signal_wind(2) signal_wind(2)],ylim,'r:');
    hold off
    legend([h2 h4])
    title(' EEEG Voltage Change whole data set')
    drawnow
    
    
    
    
    
end

end

