% script to test the accuracy of the processing code

Failed = 0;
Amp_tolerance = 1e-4; %how much error is ok for amp
Phase_tolerance = 1e-6;

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
        Amp_error_tot(iInjPhase,iPhase) = mean(Amp_error);
        Phase_error_tot(iInjPhase,iPhase) = mean(Phase_error);
        
    end
end
%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(MeasPhaseDiff,Amp_error_tot)
    ylabel('Error');
    xlabel('Phase difference beetween inj and meas');
    title('T1: Amp error for different starting phase');
    
    subplot(2,1,2)
    plot(MeasPhaseDiff,Phase_error_tot)
    ylabel('Error');
    xlabel('Phase difference');
    title('T1: Phase error for different starting phase');
    
end

Test1_AmpError=max(max(abs(Amp_error_tot)));
Test1_PhaseError=max(max(abs(Phase_error_tot)));


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
    Amp_error_tot(iAmp) = mean(Amp_error);
    Phase_error_tot(iAmp) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Amp_Meas*Amp_sc,Amp_error_tot)
    ylabel('Error');
    xlabel('Amplitude');
    title('T2: Amp error for different Amplitudes');
    
    subplot(2,1,2)
    plot(Amp_Meas*Amp_sc,Phase_error_tot)
    ylabel('Error');
    xlabel('Amplitude');
    title('T2: Phase error for different Amplitudes');
    
end

Test2_AmpError=max(max(abs(Amp_error_tot)));
Test2_PhaseError=max(max(abs(Phase_error_tot)));


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

Fc = 5:20:2000;
FreqNum = size(Fc,2);

Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_tot(iFreq) = mean(Amp_error);
    Phase_error_tot(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_tot)
    ylabel('Error');
    xlabel('Frequency');
    title('T3: Amp error 10s injection');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_tot)
    ylabel('Error');
     xlabel('Frequency');
    title('T3: Phase error 10s injection');
    
end

Test3_AmpError=max(max(abs(Amp_error_tot)));
Test3_PhaseError=max(max(abs(Phase_error_tot)));

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

Fc = 5:20:2000;
FreqNum = size(Fc,2);

Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

for iFreq = 1:FreqNum
    
    evalc('[Amp_error, Phase_error] = check_acc( Fc(iFreq),InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);');
    Amp_error_tot(iFreq) = mean(Amp_error);
    Phase_error_tot(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_tot)
    ylabel('Error');
     xlabel('Frequency');
    title('T4: Amp error 0.5s injection');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_tot)
    ylabel('Error');
     xlabel('Frequency');
    title('T4: Phase error 0.5s injection');
    
end

Test4_AmpError=max(max(abs(Amp_error_tot)));
Test4_PhaseError=max(max(abs(Phase_error_tot)));

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

Fc = 5:20:2000;
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
    Amp_error_tot(iFreq) = mean(Amp_error);
    Phase_error_tot(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_tot)
    ylabel('Error');
    xlabel('Frequency');
    title('T5: Amp error 32 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_tot)
    ylabel('Error');
    xlabel('Frequency');
    title('T5: Phase error 32 cycles');
    
end

Test5_AmpError=max(max(abs(Amp_error_tot)));
Test5_PhaseError=max(max(abs(Phase_error_tot)));

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

%% test 5

clear Amp_error_tot Phase_error_tot

%changing freq - fixed cycles of 3
disp('Test 6 - fixed cycles of 3');

Fc = 5:20:2000;
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
    Amp_error_tot(iFreq) = mean(Amp_error);
    Phase_error_tot(iFreq) = mean(Phase_error);
    
end

%%
if plotflag
    figure
    
    subplot(2,1,1)
    plot(Fc,Amp_error_tot)
    ylabel('Error');
    xlabel('Frequency');
    title('T6: Amp error 3 cycles');
    
    subplot(2,1,2)
    plot(Fc,Phase_error_tot)
    ylabel('Error');
    xlabel('Frequency');
    title('T6: Phase error 3 cycles');
    
end

Test6_AmpError=max(max(abs(Amp_error_tot)));
Test6_PhaseError=max(max(abs(Phase_error_tot)));

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











































