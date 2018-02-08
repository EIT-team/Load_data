% script to test the accuracy of the processing code

Failed = 0;
Amp_tolerance = 5e-3; %how much error is ok for amp use about 10e-5 smaller than amp
Phase_tolerance = 5e-3;

plotflag = 1;

%% test 1
% fixed amplitude, time and freq. Sweep through phase

disp('Test 1');

Fc = 1800;

InjTime=2;

Amp_Inj = 500;
Amp_Meas = 150;

InjPhase = -95:20:95;
InjPhaseNum = size(InjPhase,2);
MeasPhaseDiff = -355:40:355;
PhaseNum = size(MeasPhaseDiff,2);

Amp_error_tot=zeros(InjPhaseNum,PhaseNum);
Phase_error_tot=zeros(InjPhaseNum,PhaseNum);

for iInjPhase = 1:InjPhaseNum
    curInjPhase = InjPhase(iInjPhase);
    for iPhase = 1:PhaseNum
        
        evalc('[Amp_error, Phase_error] = check_acc( Fc,InjTime,Amp_Inj,Amp_Meas,curInjPhase,MeasPhaseDiff(iPhase) );');
        Amp_error_totT1(iInjPhase,iPhase) = mean(Amp_error);
        Phase_error_totT1(iInjPhase,iPhase) = mean(Phase_error);
        
    end
end
%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(MeasPhaseDiff,Amp_error_totT1)
    ylabel('Error');
    xlabel('Phase difference beetween inj and meas');
    title('T1: Amp error for different starting phase');
    
    subplot(2,1,2)
    plot(MeasPhaseDiff,Phase_error_totT1)
    ylabel('Error');
    xlabel('Phase difference');
    title('T1: Phase error for different starting phase');
    
end

Test1_AmpError=max(max(abs(Amp_error_totT1)));
Test1_PhaseError=max(max(abs(Phase_error_totT1)));


try
    assert( Test1_AmpError < Amp_tolerance, 'Test1 Amp error failed');
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test1_PhaseError < Phase_tolerance, 'Test1 Phase error failed');
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

%% test 2

clear Amp_error_tot Phase_error_tot

%changing amplitude
disp('Test 2');

InjTime=2;
Fc = 1800;
Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

Amp_sc=0.5:0.5:10;
AmpNum = size(Amp_sc,2);

for iAmp = 1:AmpNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc,InjTime,Amp_Inj*Amp_sc(iAmp),Amp_Meas*Amp_sc(iAmp),InjPhase,MeasPhaseDiff);');
    Amp_error_totT2(iAmp) = mean(Amp_error);
    Phase_error_totT2(iAmp) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Amp_Meas*Amp_sc,Amp_error_totT2)
    ylabel('Error');
    xlabel('Amplitude');
    title('T2: Amp error for different Amplitudes');
    
    subplot(2,1,2)
    plot(Amp_Meas*Amp_sc,Phase_error_totT2)
    ylabel('Error');
    xlabel('Amplitude');
    title('T2: Phase error for different Amplitudes');
    
end

Test2_AmpError=max(max(abs(Amp_error_totT2)));
Test2_PhaseError=max(max(abs(Phase_error_totT2)));


try
    assert( Test2_AmpError < Amp_tolerance, 'Test2 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test2_PhaseError < Phase_tolerance, 'Test2 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

%% test 3

clear Amp_error_tot Phase_error_tot

%changing freq - large fixed inj time
disp('Test 3  - large fixed injection time 10s');

InjTime=10;

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_totT3(iFreq) = mean(Amp_error);
    Phase_error_totT3(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT3)
    ylabel('Error');
    xlabel('Frequency');
    title('T3: Amp error 10s injection');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT3)
    ylabel('Error');
    xlabel('Frequency');
    title('T3: Phase error 10s injection');
    
end

Test3_AmpError=max(max(abs(Amp_error_totT3)));
Test3_PhaseError=max(max(abs(Phase_error_totT3)));

try
    assert( Test3_AmpError < Amp_tolerance, 'Test3 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test3_PhaseError < Phase_tolerance, 'Test3 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end


%% test 4

clear Amp_error_tot Phase_error_tot

%changing freq - large fixed inj time
disp('Test 4 - small fixed injection time 0.5s');

InjTime=.5;

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_totT4(iFreq) = mean(Amp_error);
    Phase_error_totT4(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT4)
    ylabel('Error');
    xlabel('Frequency');
    title('T4: Amp error 0.5s injection');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT4)
    ylabel('Error');
    xlabel('Frequency');
    title('T4: Phase error 0.5s injection');
    
end

Test4_AmpError=max(max(abs(Amp_error_totT4)));
Test4_PhaseError=max(max(abs(Phase_error_totT4)));

try
    assert( Test4_AmpError < Amp_tolerance, 'Test4 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test4_PhaseError < Phase_tolerance, 'Test4 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

%% test 5

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 32
disp('Test 5 - fixed cycles of 32');

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Cycles = 32;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_totT5(iFreq) = mean(Amp_error);
    Phase_error_totT5(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT5)
    ylabel('Error');
    xlabel('Frequency');
    title('T5: Amp error 32 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT5)
    ylabel('Error');
    xlabel('Frequency');
    title('T5: Phase error 32 cycles');
    
end

Test5_AmpError=max(max(abs(Amp_error_totT5)));
Test5_PhaseError=max(max(abs(Phase_error_totT5)));

try
    assert( Test5_AmpError < Amp_tolerance, 'Test5 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test5_PhaseError < Phase_tolerance, 'Test5 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

%% test 6

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 3
disp('Test 6 - fixed cycles of 3');

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Cycles = 3;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_totT6(iFreq) = mean(Amp_error);
    Phase_error_totT6(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT6)
    ylabel('Error');
    xlabel('Frequency');
    title('T6: Amp error 3 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT6)
    ylabel('Error');
    xlabel('Frequency');
    title('T6: Phase error 3 cycles');
    
end

Test6_AmpError=max(max(abs(Amp_error_totT6)));
Test6_PhaseError=max(max(abs(Phase_error_totT6)));

try
    assert( Test6_AmpError < Amp_tolerance, 'Test6 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test6_PhaseError < Phase_tolerance, 'Test6 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end


%% test 7

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 64
disp('Test 7 - fixed cycles of 64');

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Cycles = 64;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_totT7(iFreq) = mean(Amp_error);
    Phase_error_totT7(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT7)
    ylabel('Error');
    xlabel('Frequency');
    title('T7: Amp error 64 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT7)
    ylabel('Error');
    xlabel('Frequency');
    title('T7: Phase error 64 cycles');
    
end

Test7_AmpError=max(max(abs(Amp_error_totT7)));
Test7_PhaseError=max(max(abs(Phase_error_totT7)));

try
    assert( Test7_AmpError < Amp_tolerance, 'Test7 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test7_PhaseError < Phase_tolerance, 'Test7 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

%% test 8

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 128
disp('Test 8 - fixed cycles of 128');

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Cycles = 128;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_totT8(iFreq) = mean(Amp_error);
    Phase_error_totT8(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT8)
    ylabel('Error');
    xlabel('Frequency');
    title('T8: Amp error 128 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT8)
    ylabel('Error');
    xlabel('Frequency');
    title('T8: Phase error 128 cycles');
    
end

Test8_AmpError=max(max(abs(Amp_error_totT8)));
Test8_PhaseError=max(max(abs(Phase_error_totT8)));

try
    assert( Test8_AmpError < Amp_tolerance, 'Test8 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test8_PhaseError < Phase_tolerance, 'Test8 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end


%% test 9

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 128
disp('Test 9 - fixed cycles of 128 with dc offset');

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Cycles = 128;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 100;
DCoffsetinj = 400;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj );');
    Amp_error_totT9(iFreq) = mean(Amp_error);
    Phase_error_totT9(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT9)
    ylabel('Error');
    xlabel('Frequency');
    title('T9: Amp error 128 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT9)
    ylabel('Error');
    xlabel('Frequency');
    title('T9: Phase error 128 cycles');
    
end

Test9_AmpError=max(max(abs(Amp_error_totT9)));
Test9_PhaseError=max(max(abs(Phase_error_totT9)));

try
    assert( Test9_AmpError < Amp_tolerance, 'Test9 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test9_PhaseError < Phase_tolerance, 'Test9 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end
%% test 10

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 32
disp('Test 10 - fixed cycles of 32 with dc offset');

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Cycles = 32;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 100;
DCoffsetinj = 400;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj );');
    Amp_error_totT10(iFreq) = mean(Amp_error);
    Phase_error_totT10(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT10)
    ylabel('Error');
    xlabel('Frequency');
    title('T10: Amp error 32 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT10)
    ylabel('Error');
    xlabel('Frequency');
    title('T10: Phase error 32 cycles');
    
end

Test10_AmpError=max(max(abs(Amp_error_totT10)));
Test10_PhaseError=max(max(abs(Phase_error_totT10)));

try
    assert( Test10_AmpError < Amp_tolerance, 'Test10 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test10_PhaseError < Phase_tolerance, 'Test10 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

%% test 11

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 3
disp('Test 11 - fixed cycles of 3 with dc offset');

Fc = [5:5:150 200:200:2000];
FreqNum = size(Fc,2);

Cycles = 3;
T=(1./Fc); %Period in s
InjTime=(T.*Cycles);


Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 100;
DCoffsetinj = 400;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj );');
    Amp_error_totT11(iFreq) = mean(Amp_error);
    Phase_error_totT11(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_totT11)
    ylabel('Error');
    xlabel('Frequency');
    title('T11: Amp error 3 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_totT11)
    ylabel('Error');
    xlabel('Frequency');
    title('T11: Phase error 3 cycles');
    
end

Test11_AmpError=max(max(abs(Amp_error_totT11)));
Test11_PhaseError=max(max(abs(Phase_error_totT11)));

try
    assert( Test11_AmpError < Amp_tolerance, 'Test11 Amp error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

try
    assert( Test11_PhaseError < Phase_tolerance, 'Test11 Phase error failed')
catch err
    Failed = 1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

% SNR TEST


%%

figure
hold on
plot(Fc,Amp_error_totT5);
plot(Fc,Amp_error_totT7);
plot(Fc,Amp_error_totT8);
hold off
ylabel('Error');
xlabel('Frequency');
legend('32','64','128')




%% Check a given ExpSetup

load('S3b_MF1_log.mat','ExpSetup');



clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 3
disp('Test 12 - ExpSetup');

Fc = ExpSetup.Freq';
FreqNum = size(Fc,2);


InjTime=ExpSetup.MeasurementTime/1000';


Amp_Inj = 9500;
Amp_Meas = 1200 ;
InjPhase=0;
MeasPhaseDiff=-30;
DCoffset = 100;
DCoffsetinj = 400;

for iFreq = 1:FreqNum
    
    %     evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj );');
    [Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime(iFreq),Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff,DCoffset,DCoffsetinj );
    
    
    Amp_error_totT12(iFreq) = mean(Amp_error);
    Phase_error_totT12(iFreq) = mean(Phase_error);
    
end

%%

if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,100*(Amp_error_totT12/Amp_Inj))
    ylabel('Error %');
    xlabel('Frequency');
    title('T12: Amp error 3 cycles');
    
    subplot(2,1,2)
    
    ylabel('Error');
    xlabel('Frequency');
    title('T12: Phase error 3 cycles');
    
end
%%

figure
hold on
plot(Fc/1000,100*(Amp_error_totT12/Amp_Inj))
plot(Fc/1000,100*(Phase_error_totT12/MeasPhaseDiff))

ylabel('Error %');
xlabel('Frequency (kHz)');

legend('Amplitude','Phase');



%%




















