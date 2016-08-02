%%


Fc=5;

Fpass =2;
Fstop = 0.2;
Apass = 0.5;
Astop =10;
Fs = 16384;

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
    implen=impzlength(d,Decay_coef);
    
    
    [H,W] =freqz(d,Fc-10:0.2:Fc+10,Fs);
fc_idx = find (W > Fc,1);
GainFc =abs(H(fc_idx));

ok = 1- GainFc < tol;