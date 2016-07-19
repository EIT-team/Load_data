Fc = 500;

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
BW=100;
StopBandDiff=400;

Fc = 500;

Fpass1 = Fc-BW/2;
Fstop1 = Fpass1 -StopBandDiff;

Fpass2 = Fc+BW/2;
Fstop2 = Fpass2 + StopBandDiff;


Astop1 = 60;
Astop2= Astop1;
Apass  = 0.5;


Wp = [Fpass1 Fpass2]/Fs/2;
Ws = [Fstop1 Fstop2]/Fs/2;
Rp = Apass;
Rs = Astop1;
[n,Wn] = buttord(Wp,Ws,Rp,Rs)


h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
    Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);

d = design(h, 'butter', ...
    'MatchExactly', 'passband', ...
    'SOSScaleNorm', 'Linf');

st= isstable(d);
imp_len=impzlength(d,0.001);
imp_len_ms = (imp_len / Fs)*(1000);
len_cyc = (1/Fc)*1000;
imp_len_cyc=imp_len_ms/len_cyc;

disp([' istable : ' num2str(st)]);
disp([' impzlength : ' num2str(imp_len)]);
disp([' impzlength ms: ' num2str(imp_len_ms)]);
disp([' impzlength cycl: ' num2str(imp_len_cyc)]);

[B,A] = butter(3,(Fc+[-BW/2,BW/2])./(Fs/2));
% [B,A] = butter(n,Wn);

st= isstable(B,A);
imp_len=impzlength(B,A,0.001);
imp_len_ms = (imp_len / Fs)*(1000);
len_cyc = (1/Fc)*1000;
imp_len_cyc=imp_len_ms/len_cyc;

disp([' istable : ' num2str(st)]);
disp([' impzlength : ' num2str(imp_len)]);
disp([' impzlength ms: ' num2str(imp_len_ms)]);
disp([' impzlength cycl: ' num2str(imp_len_cyc)]);



[H1,W]  = freqz(d,Fstop1:0.2:Fstop2,Fs);
[H2,W] =freqz(B,A,Fstop1:0.2:Fstop2,Fs);

figure;
hold on
plot(W,10*log10(abs(H1)))
plot(W,10*log10(abs(H2)))
hold off
legend('New','Old')


%%

Fs=16384;
BW=15;
StopBandDiff=150;

Fc = 25;

%high pass
BWhp = 2;
StopBandDiffhp=15;
Fpass1 = Fc-BWhp/2;
Fstop1 = Fpass1 -StopBandDiffhp;

%low pass
Fpass2 = Fc+BW/2;
Fstop2 = Fpass2 + StopBandDiff;


Astop1 = 60;
Astop2= Astop1;
Apass  = 0.5;


% hhp = fdesign.highpass('fp,fst,ap,ast', Fpass1, Fstop1, Apass, Astop1, Fs);
% 
% dhp = design(hhp, 'butter', ...
%     'MatchExactly', 'stopband', ...
%     'SOSScaleNorm', 'Linf');

hhp = fdesign.highpass('fst,fp,ast,ap', Fstop1, Fpass1, Astop1, Apass, Fs);

dhp = design(h, 'butter', ...
    'MatchExactly', 'passband', ...
    'SOSScaleNorm', 'Linf');

st= isstable(dhp);
imp_len=impzlength(dhp,0.001);
imp_len_ms = (imp_len / Fs)*(1000);
len_cyc = (1/Fc)*1000;
imp_len_cyc=imp_len_ms/len_cyc;

disp(['hp istable : ' num2str(st)]);
disp(['hp impzlength : ' num2str(imp_len)]);
disp(['hp impzlength ms: ' num2str(imp_len_ms)]);
disp(['hp impzlength cycl: ' num2str(imp_len_cyc)]);


hlp = fdesign.lowpass('fp,fst,ap,ast', Fpass2, Fstop2, Apass, Astop2, Fs);

dlp = design(hlp, 'butter', ...
    'MatchExactly', 'passband', ...
    'SOSScaleNorm', 'Linf');

st= isstable(dlp);
imp_len=impzlength(dlp,0.001);
imp_len_ms = (imp_len / Fs)*(1000);
len_cyc = (1/Fc)*1000;
imp_len_cyc=imp_len_ms/len_cyc;

disp(['lp istable : ' num2str(st)]);
disp(['lp impzlength : ' num2str(imp_len)]);
disp(['lp impzlength ms: ' num2str(imp_len_ms)]);
disp(['lp impzlength cycl: ' num2str(imp_len_cyc)]);




[B,A] = butter(3,(Fc+[-BW/2,BW/2])./(Fs/2));
% [B,A] = butter(n,Wn);

st= isstable(B,A);
imp_len=impzlength(B,A,0.001);
imp_len_ms = (imp_len / Fs)*(1000);
len_cyc = (1/Fc)*1000;
imp_len_cyc=imp_len_ms/len_cyc;

disp([' old istable : ' num2str(st)]);
disp([' old impzlength : ' num2str(imp_len)]);
disp([' old impzlength ms: ' num2str(imp_len_ms)]);
disp([' old impzlength cycl: ' num2str(imp_len_cyc)]);



[H1,W]  = freqz(dlp,Fstop1:0.2:Fstop2,Fs);
[H2,W] =freqz(B,A,Fstop1:0.2:Fstop2,Fs);
[H3,W] =freqz(dhp,Fstop1:0.2:Fstop2,Fs);

figure;
hold on
plot(W,10*log10(abs(H1)))
plot(W,10*log10(abs(H2)))
plot(W,10*log10(abs(H3)))
hold off
legend('New','Old','newhp')



