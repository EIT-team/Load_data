Fcur = 15;
FreqNum = size(Fcur,2);

Cycles = 128;
T=(1./Fcur); %Period in s
InjTime=(T.*Cycles);

% InjTime=10;


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 0;
DCoffsetinj = 0;


[Amp_error1, Phase_error1,V1,Vd1,Filt1,tr1] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj,[]);

trim_demod=200;
Fs=16384;

%%

Fcur = 15;
FreqNum = size(Fcur,2);

Cycles = 128;
T=(1./Fcur); %Period in s
InjTime=(T.*Cycles);

% InjTime=10;


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 100;
DCoffsetinj = 400;


[Amp_error2, Phase_error2,V2,Vd2,Filt2,tr2] = check_acc( Fcur,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj,[]);



fvtool(Filt1,Filt2)


%

%%
figure
hold on
plot(V1)
plot(V2)
% plot(V8)
hold off


figure
hold on
% plot(Vd1(tr1:end-tr1))
plot(Vd1)
plot(Vd2)
% plot(Vd2(tr2:end-tr2))
% plot(Vd8)
hold off
% ylim(Amp_Inj + [-1 +1])


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