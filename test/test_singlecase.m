Fc = 25;

Cycles = 3;
T=(1./Fc); %Period in s
InjTime=ceil(T.*Cycles);


% Amp_Inj = 500;
% Amp_Meas = 150;
% InjPhase=0;
% MeasPhaseDiff=-30;
% 
% 
% [Amp_error, Phase_error] = check_acc( Fc,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff);


%%

Fs=16384;
BW=50;
StopBandDiff=250;

Fc = 500;

[B,A] = butter(3,(Fc+[-BW/2,BW/2])./(Fs/2));

Fpass1 = Fc-BW/2;
Fstop1 = Fpass1 -StopBandDiff;

Fpass2 = Fc+BW/2;
Fstop2 = Fpass2 + StopBandDiff;


Astop1 = 60;
Apass  = 0.5;
Astop2 = 60;

% d = designfilt('bandpassiir', ...
%   'StopbandFrequency1',Fstop1,'PassbandFrequency1', Fpass1, ...
%   'PassbandFrequency2',Fpass2,'StopbandFrequency2', Fstop2, ...
%   'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
%   'StopbandAttenuation2',Astop2, ...
%   'DesignMethod','butter','SampleRate', Fs);



h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
    Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);

d = design(h, 'butter', ...
    'MatchExactly', 'passband', ...
    'SOSScaleNorm', 'Linf');


 fvtool(d)

[b,a]=d.tf;

st= isstable(d);
imp_len=impzlength(d,0.001);
imp_len_ms = (imp_len / Fs)*(1000);

disp([' istable : ' num2str(isstable(d))]);
disp([' impzlength : ' num2str(imp_len)]);
disp([' impzlength ms: ' num2str(imp_len_ms)]);


[H1,W]  = freqz(d,10000,Fs);
[H2,W] =freqz(B,A,10000,Fs);

figure;
hold on
plot(W,10*log10(abs(H1)))
plot(W,10*log10(abs(H2)))
hold off
legend('New','Old')

%% low freq ones - two filers?

% Low pass one
Fpass = 50;     % Passband Frequency
Fstop = 200;    % Stopband Frequency
Apass = 0.5;    % Passband Ripple (dB)
Astop = 60;     % Stopband Attenuation (dB)
Fs    = 16384;  % Sampling Frequency

h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);

Hd = design(h, 'butter', ...
    'MatchExactly', 'stopband', ...
    'SOSScaleNorm', 'Linf');

%%

Fs=16384;
BW=5;
StopBandDiff=10;

Fc = 25;

Fpass1 = Fc-BW/2;
Fstop1 = Fpass1 -StopBandDiff;

Fpass2 = Fc+BW/2;
Fstop2 = Fpass2 + StopBandDiff;


Astop1 = 60;
Apass  = 0.5;


Wp = [Fpass1 ]/Fs/2;
Ws = [Fstop1 ]/Fs/2;
Rp = Apass;
Rs = Astop1;
[n,Wn] = buttord(Wp,Ws,Rp,Rs)









