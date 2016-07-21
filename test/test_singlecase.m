Fc = 55;
FreqNum = size(Fc,2);

Cycles = 32;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);

Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

[Amp_error, Phase_error] = check_acc( Fc,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);


%


