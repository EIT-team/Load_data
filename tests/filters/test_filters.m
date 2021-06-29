%%

% Fs=16384;
Fs=100000;

Fc_all = 5:100:550;
FreqNum = length(Fc_all);

BWtarget =100;

Fstopomin=0.5;
Fcomin = 1.5;


for iFreq = 1: FreqNum
    
    
    Fc=Fc_all(iFreq);
    
    disp(Fc);
    
    if Fc == 25;
        disp('stop');
    end
    
    %% old filter
    
    BWold = BWtarget;
    Fco(2) = Fc+(BWold/2);
    Fco(1) = max([(Fc -(BWold/2)), Fstopomin]);
    [B,A] = butter(3,Fco./(Fs/2));
    
    stold(iFreq)= isstable(B,A);
    imp_lenold(iFreq)=impzlength(B,A,0.001);
    imp_len_msold(iFreq) = (imp_lenold(iFreq) / Fs)*(1000);
    len_cyc = (1/Fc)*1000;
    imp_len_cycold(iFreq)=imp_len_msold(iFreq)/len_cyc;
    
    [Hold,W] =freqz(B,A,Fco(1):0.2:Fco(2),Fs);
    fc_idx = find (W > Fc,1);
    Gainold(iFreq) =abs(Hold(fc_idx));
    
    
    %% New Bandpass
    
    BWbp = BWtarget;
    
    StopBandDiffbp = 350;
    
    Fpass1bp = Fc-BWbp/2;
    
    Astop1 = 60;
    Astop2= 60;
    Apass  = 0.5;
    
    
    if ~(Fpass1bp > 0)
        Fpass1bp = max([Fc/2 Fcomin]);
        Astop1 = 15;
    end
    
    Fstop1bp = Fpass1bp -StopBandDiffbp;
    Fstop1bp = max([Fstop1bp Fstopomin]);
    
    Fpass2bp = Fc+BWbp/2;
    Fstop2bp = Fpass2bp + StopBandDiffbp;
    
    
    
    hbp = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1bp, Fpass1bp, ...
        Fpass2bp, Fstop2bp, Astop1, Apass, Astop2, Fs);
    
    dbp = design(hbp, 'butter', 'MatchExactly', 'passband', 'SOSScaleNorm', 'Linf');
    
    %     dbp = designfilt('bandpassiir', ...       % Response type
    %         'FilterOrder', 20, ...
    %        'HalfPowerFrequency1',Fpass1bp, ...    % Frequency constraints
    %        'HalfPowerFrequency2',Fpass2bp, ...
    %        'SampleRate',Fs)   ;            % Sample rate
    %
    
    
    
    stbp(iFreq)= isstable(dbp);
    imp_lenbp(iFreq)=impzlength(dbp,0.001);
    imp_len_msbp(iFreq) = (imp_lenbp(iFreq) / Fs)*(1000);
    len_cyc = (1/Fc)*1000;
    imp_len_cycbp(iFreq)=imp_len_msbp(iFreq)/len_cyc;
    
    [Hbp,W] =freqz(dbp,Fstop1bp:0.2:Fstop2bp,Fs);
    fc_idx = find (W > Fc,1);
    Gainbp(iFreq) =abs(Hbp(fc_idx));
    
    
    %% New Cascaded filter
    
    BWc = BWtarget;
    
    StopBandDiffc = 350;
    
    Fpass1c = Fc-BWc/2;
    
    Astop1 = 60;
    Astop2= 60;
    Apass  = 0.5;
    
    if ~(Fpass1c > 0)
        Fpass1c = max([Fc/2 Fcomin]);
        Astop1 =30;
    end
    
    Fstop1c = Fpass1c -StopBandDiffc;
    Fstop1c = max([Fstop1c Fstopomin]);
    
    Fpass2c = Fc+BWc/2;
    Fstop2c = Fpass2c + StopBandDiffc;
    
    
    
    %     hchp = fdesign.highpass('fst,fp,ast,ap', Fstop1c, Fpass1c, Astop1, Apass, Fs);
    %
    %     dchp = design(hchp, 'butter', 'MatchExactly', 'passband', 'SOSScaleNorm', 'Linf');
    
    %     hclp = fdesign.lowpass('fp,fst,ap,ast', Fpass2c, Fstop2c, Apass, Astop2, Fs);
    %
    %     dclp = design(hclp, 'butter', 'MatchExactly', 'passband', 'SOSScaleNorm', 'Linf');
    
    dclp = designfilt('lowpassiir', ...       % Response type
        'FilterOrder', 10, ...
        'HalfPowerFrequency',Fpass2c, ...    % Frequency constraints
        'SampleRate',Fs)   ;            % Sample rate
    
    
    dchp = designfilt('highpassiir', ...       % Response type
        'FilterOrder', 10, ...
        'HalfPowerFrequency',Fpass1c, ...    % Frequency constraints
        'SampleRate',Fs)   ;            % Sample rate
    
    
    stc(iFreq)= all([isstable(dchp) isstable(dclp)]);
    imp_lenclp(iFreq)=impzlength(dclp,0.001) ;
    imp_lenchp(iFreq) = impzlength(dchp,0.001);
    imp_len_msclp(iFreq) = (imp_lenclp(iFreq) / Fs)*(1000);
    imp_len_mschp(iFreq) = (imp_lenchp(iFreq) / Fs)*(1000);
    len_cyc = (1/Fc)*1000;
    imp_len_cycclp(iFreq)=imp_len_msclp(iFreq)/len_cyc;
    imp_len_cycchp(iFreq)=imp_len_mschp(iFreq)/len_cyc;
    [Hclp,W] =freqz(dclp,Fstop1c:0.2:Fstop2c,Fs);
    [Hchp,W] =freqz(dchp,Fstop1c:0.2:Fstop2c,Fs);
    fc_idx = find (W > Fc,1);
    Gainclp(iFreq) = abs(Hclp(fc_idx)) ;
    
    Gainchp(iFreq) =  abs(Hchp(fc_idx));
    
    %% Freq Info
    
end


%%

figure;
hold on
plot(Fc_all,10*log10(Gainold));
plot(Fc_all,10*log10(Gainbp));
plot(Fc_all,10*log10(Gainclp));
plot(Fc_all,10*log10(Gainchp));
hold off
legend('Old','Bandp','Casclp','Caschp')
% legend('Bandp','Casclp','Caschp')
title('gain')
xlabel('Freq')
ylabel('Gain db');


figure;
hold on
plot(Fc_all,imp_len_msold);
plot(Fc_all,imp_len_msbp);
plot(Fc_all,imp_len_msclp);
plot(Fc_all,imp_len_mschp);

hold off
legend('Old','Bandp','Casclp','Caschp')
% legend('Bandp','Casclp','Caschp')
title('timedelay ms')
ylabel('Settling time ms')
xlabel('Freq')



figure;
hold on
plot(Fc_all,imp_len_cycold);
plot(Fc_all,imp_len_cycbp);
plot(Fc_all,imp_len_cycclp);
plot(Fc_all,imp_len_cycchp);

hold off
legend('Old','Bandp','Casclp','Caschp')
% legend('Bandp','Casclp','Caschp')
title('timedelay cycles')
ylabel('Settling time cycles')
xlabel('Freq')





