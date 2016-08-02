Fcur = 15;
FreqNum = size(Fcur,2);

Cycles = 64;
T=(1./Fcur); %Period in s
InjTime=(T.*Cycles);

% InjTime=10;


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 0;
DCoffsetinj = 0;

Fs=100000;


[Amp_error1, Phase_error1,V1,Vd1,Filt1,tr1] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj,[],Fs);

trim_demod=200;


t1 = (0:1:length(V1)-1)/Fs;


%%

Fcur = 20;
FreqNum = size(Fcur,2);

Cycles = 64;
T=(1./Fcur); %Period in s
InjTime=(T.*Cycles);

% InjTime=10;


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 0;
DCoffsetinj = 0;

Fs2=100000;


[Amp_error2, Phase_error2,V2,Vd2,Filt2,tr2] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj,[],Fs2);


t2 = (0:1:length(V2)-1)/Fs2;

%

%%

fvtool(Filt1,Filt2)
%%
[H1,F1]=freqz(Filt1,Fcur-500:Fcur+500,Fs);
[H2,F2]=freqz(Filt2,Fcur-500:Fcur+500,Fs2);

figure;
hold on
plot(F1,10*log10(abs(H1)));
plot(F2,10*log10(abs(H2)));
hold off

%%

figure
hold on
plot(t1,V1)
plot(t2,V2)
% plot(V8)
hold off


figure
hold on
% plot(Vd1(tr1:end-tr1))
plot(t1,Vd1)
plot(t2,Vd2)
% plot(Vd2(tr2:end-tr2))
% plot(Vd8)
hold off
% ylim(Amp_Inj + [-1 +1])
