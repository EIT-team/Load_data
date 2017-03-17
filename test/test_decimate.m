Fcur = 1000;
FreqNum = size(Fcur,2);

% Cycles = 2*6000;
% T=(1./Fcur); %Period in s
% InjTime=(T.*Cycles);

InjTime=2;


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=0;
DCoffset = 0;
DCoffsetinj = 0;

Fs=100000;


[Amp_error1, Phase_error1,V1,Vd1,Filt1,tr1] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj,[],Fs);
%%
dec_factor=2;

Fs2=Fs/dec_factor;

[Amp_error2, Phase_error2,V2,Vd2,Filt2,tr2] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj,[],Fs2);
%%
[Amp_error3, Phase_error3,V3,Vd3,Filt3,tr3] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj,[],Fs,dec_factor);



%%

t1 = (0:1:length(Vd1)-1)/Fs;
t2 = (0:1:length(Vd2)-1)/Fs2;
t3 = (0:1:length(Vd3)-1)/Fs2;

figure;
hold on
plot(t1,Vd1);
plot(t2,Vd2);
plot(t3,Vd3);
hold off
legend('Original','Ideal downsample','Decimated');

figure;
hold on
plot(t1,V1);
plot(t2,V2);
plot(t3,V3);
hold off
legend('Original','Ideal downsample','Decimated');
xlim([1, 1.01])

