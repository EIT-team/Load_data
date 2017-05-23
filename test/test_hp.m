%%
BVsold=load('E:\testperchn\bignir\MF_run1-BVOLD');
B=BVsold.info.B{1}{2};
A=BVsold.info.A{1}{2};



%%
Fc=5;

Fpass =Fc -2.5;
Fstop = Fpass -1.5;
Apass = 0.5;
Astop =5;
Fs = 16384;
% Fs = 100000;

d = designfilt('highpassiir', ...       % Response type
    'PassbandFrequency',Fpass, ...
    'StopbandFrequency',Fstop, ...
    'StopbandAttenuation',Astop, ...   % Magnitude constraints
    'PassbandRipple',Apass, ...
    'DesignMethod','butter', ...      % Design method
    'MatchExactly','passband', ...   % Design method options
    'SampleRate',Fs)     ;

%%

st = isstable(d);
implen =impzlength(d);
implenold =impzlength(B,A);

len_cyc = (1/Fc)*Fs;
numcyc =implen/len_cyc

fprintf('imp len old %d, new %d',implenold,implen);


tol = 1e-2;

[Hold,Wold] =freqz(B,A,0:.1:Fc+10,Fs);
[Hnew,Wnew] =freqz(d,0:.1:Fc+10,Fs);

figure
hold on
    plot(Wnew,10*log10(abs(Hnew)));
    plot(Wold,10*log10(abs(Hold)));
hold off
 legend('new','old');
fc_idx = find (Wnew > Fc,1);
GainFc =abs(Hnew(fc_idx))

ok = 1- GainFc < tol