Fcur = 50;
FreqNum = size(Fcur,2);

Cycles = 32;
T=(1./Fcur); %Period in s
InjTime=(T.*Cycles);

InjTime=10;


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

[Amp_error1, Phase_error1,V1,Vd1,Filt1] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,[],[],1);

trim_demod=200;
Fs=16384;

%%

Fcur = 55;
FreqNum = size(Fcur,2);

Cycles = 32;
T=(1./Fcur); %Period in s
InjTime=(T.*Cycles);

InjTime=10;


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

[Amp_error2, Phase_error2,V2,Vd2,Filt2] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,[],[],1);



fvtool(Filt1,Filt2)


%

%%
figure
hold on
plot(V1)
% plot(V4)
% plot(V8)
hold off


figure
hold on
plot(Vd1)
% plot(Vd4)
% plot(Vd8)
hold off

%
% lpFilt = designfilt('lowpassfir','PassbandFrequency',0.25, ...
%          'StopbandFrequency',0.35,'PassbandRipple',0.5, ...
%          'StopbandAttenuation',65,'DesignMethod','kaiserwin');


%      lpFilt = designfilt('lowpassiir', ...        % Response type
%        'PassbandFrequency',50, ...     % Frequency constraints
%        'StopbandFrequency',100, ...
%        'DesignMethod','butter', ...         % Design method
%        'PassbandRipple',0.5, ...          % Design method options
%        'StopbandAttenuation',65, ...
%        'SampleRate',Fs) ;              % Sample rate


%%
[B,A]= fir1(100,5/(Fs/2));


Vs=abs(hilbert(V1));

figure;
hold on
plot(Vd1)

plot(filtfilt(B,A,Vd1));
hold off




%%
% 
% 
% HDR=ScouseTom_getHDR('E:\testperchn\bignir\MF_run1.bdf');
% Trigger= ScouseTom_TrigReadChn(HDR);
% TT=ScouseTom_TrigProcess(Trigger,HDR);
% Veeg=sread(HDR,15,0);
% %%
% vsig=Veeg(TT.InjectionSwitches{1}(1,1):TT.InjectionSwitches{1}(1,2),:);
% [ trim_demod,FilterOut,Fc ] = ScouseTom_data_GetFilterTrim( vsig,Fs);
% 
% 
% [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( vsig,FilterOut);