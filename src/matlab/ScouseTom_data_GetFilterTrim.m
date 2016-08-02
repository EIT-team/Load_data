function [ trim_demod,FilterOut,Fc] = ScouseTom_data_GetFilterTrim( Vseg,Fs,BWtarget,MaxImpSamples,Fc,plotflag )
%ScouseTom_GetFilterTrim Gets optimal filter parameters and samples to trim
%   based on SNR tests for different windows
% from the dexterous exploratory hands of jimmy




if  exist('Fc','var') == 0 || isempty(Fc)
    %get carrier frequency
    Fc=ScouseTom_data_GetCarrier(Vseg,Fs);
end



%get number of samples in segment
Nsamples=length(Vseg);


if  exist('BWtarget','var') == 0 || isempty(BWtarget)
    BWtarget =100;
end

if  exist('MaxImpSamples','var') ==0 || isempty(MaxImpSamples)
    %trim max is 25% of signal - to give 50% for the average
    MaxImpSamples=ceil(0.4*Nsamples);
    %round up to nearest 10 samples - this is so that data with window size 1
    %sample different (as can happen) produces the same filter, and thus the
    %same result.
    rndnum=10;
    MaxImpSamples=round(MaxImpSamples/rndnum)*rndnum;
    
    %dont use this corrected version if the data is very small
    if MaxImpSamples > ceil(0.5*Nsamples)
        
        MaxImpSamples=ceil(0.4*Nsamples);
        
    end
    
end


if  exist('plotflag','var') == 0 || isempty(plotflag)
    plotflag =0;
end

%amount filter ripple must decay by before using data%this is based on
%number of samples in chapter 3 in y thesis. although it could be more
%rigourously chosen
Decay_coef=0.0001;


Fstopomin=0.5;
FcominIIR = 15;
Astop = 60;
StopBandDiffMax = 350;
StopBandDiffMin = 50;

FcominFIR =8;


% Low pass only cutoff
FIRLowPassCutoff = 20;
IIRLowPassCutoff = 16;
UseLowPassIIR=0;
UseLowPassFIR=0;


%the filter takes some time in seconds to decay to ~zero for IIR butterworth filters.
%This is independent of sampling rate, so we need to have a different max
%number of samples for different sampling rates. This is the max for FIR

Mindecay_seconds=0.2;
MaxFirOrder = floor(Mindecay_seconds*Fs);


%% Find best IIR filter
% if IIR decay time this is too long then we have to use much slow FIR
% filters to avoid artefacts

IIRfailed =0;

dprev=[];

try
    
    BWIIR = BWtarget;
    
    IIRfound =0;
    
    curStopBandDiff = StopBandDiffMax;
    
    while ~IIRfound
        
        
        if Fc > IIRLowPassCutoff
            [dcur,st,gainok,implen]=BandPassIIR(BWIIR,curStopBandDiff,Astop,Fc,FcominIIR,Fstopomin,Fs,Decay_coef);
        else
            UseLowPassIIR=1;
            [dcur,st,gainok,implen]=LowPassIIR(BWIIR,curStopBandDiff,Astop,Fc,Fs,Decay_coef);
        end
        
        if st && implen < MaxImpSamples
            
            curStopBandDiff = curStopBandDiff - 50;
            dprev = dcur;
            
            if curStopBandDiff < StopBandDiffMin || ~gainok
                IIRfound =1;
            end
            
        else
            
            IIRfound =1;
            
        end
        
        
    end
    
    dIIR = dprev;
    
    if isempty(dIIR)
        IIRfailed =1;
        Samples_needed = inf;
        IIRgainok =0;
    else
        
        Samples_needed = impzlength(dIIR,Decay_coef);
        
    end
    
    
    if plotflag ==1 && ~IIRfailed;
        figure;
        hold on
        
        [H,T]=impz(dIIR,max([MaxImpSamples, implen]));
        
        plot(T,H)
        %line([0 length(T)],[maxh*Decay_coef maxh*Decay_coef],'color','r');
        line([Samples_needed Samples_needed],[min(H) max(H)],'color','r');
        line([MaxImpSamples MaxImpSamples],[min(H) max(H)],'color',[0 0.5 0]);
        hold off
        title('impulse response of IIR filter')
        %         xlim([0 (ceil(Samples_needed/1000))*1000])
        % set(gca,'Yscale','log');
        legend('Filter response','Req. Trim Samples','Max Trim Samples')
        drawnow
    end
    
catch err
    fprintf(2, '%s\n', getReport(err, 'extended'));
    IIRfailed =1 ;
end

%% Find best FIR filter

BWFIR = BWtarget;

N = min([MaxImpSamples MaxFirOrder]);

if Fc > FIRLowPassCutoff
    [dFIR,FIRgainok]=BandPassFIR(BWFIR,N,Fc,FcominFIR,Fs);
else
    [dFIR,FIRgainok]=LowPassFIR(BWFIR,N,Fc,Fs);
    UseLowPassFIR=1;
end


if Fc < 1
    ForceFIR =1;
else
    ForceFIR =0;
end


%% choose filter

%from TestFilterSNR - FIR outperforms IIR until around 1200 samples (on biosemi Fs and BW50), then
%they are within 1% of each other. IIR is MUCH faster than high order FIR
%so use this to speed up the process

if (MaxImpSamples <Samples_needed) || IIRfailed || ForceFIR
    %if we do not have enough samples to allow for the filter to decay
    %sufficiently, then use the slower FIR filter. Blackman harris window
    %chosen as it gives the best trade off between stopband ripple and
    %rolloff
    
    FilterOut = dFIR;
    
    
    trim_demod=N;
    
    if UseLowPassFIR
        
        disp('Low Pass FIR with Blackman-Harris Window used');
    else
        disp('FIR with Blackman-Harris Window used');
        
    end
else
    
    FilterOut = dIIR;
    
    %if we have more samples than we need, still use the max as we want the
    %filter to decay as much as possible. This reduces errors from ripples
    %caused by hilbert
    
    trim_demod=MaxImpSamples;
    
    
    if UseLowPassIIR
        disp('Low Pass Min Order Butterworth Filter Used');
    else
        disp('Min Order Butterworth Filter Used');
    end
    
    
end
%%
if plotflag ==1;
    %     figure;
    impz(FilterOut);
    drawnow
    %     figure;
    freqz(FilterOut,2000,Fs);
    drawnow
end

%% final gain check

[gainok,GainFc]=CheckGainFc(Fc,Fs,FilterOut,1e-3);

if ~gainok
    fprintf(2,'WARNING! GAIN NOT OK AT CARRIER FREQ! Gain : %.6f\n',GainFc);
end

if UseLowPassFIR || UseLowPassIIR
    disp('High Pass Min Order Butterworth Filter Used');
    
    FilterOut={FilterOut};
    FilterOut{2}=HighPassIIR(Fc,Fs);
    
end



end


function [d,st,gainok,implen]=BandPassIIR(BWin,StopBandDiff,Astop,Fc,Fcomin,Fstopomin,Fs,Decay_coef)


Apass  = 0.5;

gainok =0;
curBW = BWin;
iterations =0;
BWincrement = 25;
GainFcprev =0;

maxiterations =10;

while gainok ==0 && iterations < maxiterations
    
    Fpass1bp = Fc-curBW/2;
    
    Astop1 = Astop;
    Astop2= Astop;
    
    if (Fpass1bp <= Fcomin)
        Fpass1bp = max([Fc/2 Fcomin]);
        Astop1 = 15;
    end
    
    Fstop1bp = Fpass1bp -StopBandDiff;
    Fstop1bp = max([Fstop1bp Fstopomin]);
    
    Fpass2bp = Fc+curBW/2;
    Fstop2bp = Fpass2bp + StopBandDiff;
    
    %     hbp = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1bp, Fpass1bp, ...
    %         Fpass2bp, Fstop2bp, Astop1, Apass, Astop2, Fs);
    
    %     d = design(hbp, 'butter', 'MatchExactly', 'passband', 'SOSScaleNorm', 'Linf');
    
    
    d = designfilt('bandpassiir', ...       % Response type
        'StopbandFrequency1',Fstop1bp, ...    % Frequency constraints
        'PassbandFrequency1',Fpass1bp, ...
        'PassbandFrequency2',Fpass2bp, ...
        'StopbandFrequency2',Fstop2bp, ...
        'StopbandAttenuation1',Astop1, ...   % Magnitude constraints
        'PassbandRipple',Apass, ...
        'StopbandAttenuation2',Astop2, ...
        'DesignMethod','butter', ...      % Design method
        'MatchExactly','passband', ...   % Design method options
        'SampleRate',Fs)     ;          % Sample rate
    
    
    st = isstable(d);
    implen=impzlength(d,Decay_coef);
    
    [gainok,GainFc]=CheckGainFc(Fc,Fs,d,1e-2);
    
    if ~gainok
        curBW = curBW + BWincrement;
        if GainFcprev > GainFc
            d = dprev;
            gainok=0;
            fprintf(2,'warning, gain going wrong way!\n');
            iterations = maxiterations+1;
            
        else
            GainFcprev=GainFc;
            dprev =d;
        end
    end
    iterations =iterations +1;
    
end

end

function [d,st,gainok,implen]=LowPassIIR(BWin,StopBandDiff,Astop,Fc,Fs,Decay_coef)


Apass  = 0.5;

gainok =0;
curBW = BWin;
iterations =0;
BWincrement = 25;
GainFcprev =0;
while gainok ==0 && iterations < 10
    
    Fpass = Fc+curBW/2;
    Fstop= Fpass + StopBandDiff;
    
    %     hbp = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1bp, Fpass1bp, ...
    %         Fpass2bp, Fstop2bp, Astop1, Apass, Astop2, Fs);
    
    %     d = design(hbp, 'butter', 'MatchExactly', 'passband', 'SOSScaleNorm', 'Linf');
    
    
    d = designfilt('lowpassiir', ...       % Response type
        'PassbandFrequency',Fpass, ...
        'StopbandFrequency',Fstop, ...
        'StopbandAttenuation',Astop, ...   % Magnitude constraints
        'PassbandRipple',Apass, ...
        'DesignMethod','butter', ...      % Design method
        'MatchExactly','passband', ...   % Design method options
        'SampleRate',Fs)     ;          % Sample rate
    
    
    st = isstable(d);
    implen=impzlength(d,Decay_coef);
    
    [gainok,GainFc]=CheckGainFc(Fc,Fs,d,1e-2);
    
    if ~gainok
        curBW = curBW + BWincrement;
        if GainFcprev > GainFc
            d = dprev;
            gainok=1;
        else
            GainFcprev=GainFc;
            dprev =d;
        end
    end
    iterations =iterations +1;
    
end

end

function [d]=HighPassIIR(Fc,Fs)

Fpass = max([(Fc -2.5) 1]);
Fstop = max([(Fpass -1.5) 0.5]);
Apass = 0.5;
Astop =5;

d = designfilt('highpassiir', ...       % Response type
    'PassbandFrequency',Fpass, ...
    'StopbandFrequency',Fstop, ...
    'StopbandAttenuation',Astop, ...   % Magnitude constraints
    'PassbandRipple',Apass, ...
    'DesignMethod','butter', ...      % Design method
    'MatchExactly','passband', ...   % Design method options
    'SampleRate',Fs)     ;

end




function [d,gainok]=BandPassFIR(BWin,N,Fc,FcominFIR,Fs)

gainok =0;
curBW = BWin;
iterations =0;
BWincrement = 25;
GainFcprev =0;
while gainok ==0 && iterations < 10
    
    F6dB1 = Fc-curBW/2;
    
    
    if ~(F6dB1 > FcominFIR)
        %         F6dB1 = max([Fc/2 Fcomin]);
        F6dB1 = max([FcominFIR]);
    end
    
    
    F6dB2 = Fc+curBW/2;
    
    
    %     h = fdesign.bandpass('n,fc1,fc2', N, F6dB1, F6dB2, Fs);
    %
    %     d = design(h, 'window', 'Window', 'blackmanharris');
    
    d = designfilt('bandpassfir', ...       % Response type
        'FilterOrder',N, ...            % Filter order
        'CutoffFrequency1',F6dB1, ...    % Frequency constraints
        'CutoffFrequency2',F6dB2, ...
        'DesignMethod','window', ...         % Design method
        'Window','blackmanharris', ...         % Design method options
        'SampleRate',Fs);               % Sample rate
    
    
    
    [gainok,GainFc]=CheckGainFc(Fc,Fs,d,1e-2);
    
    if ~gainok
        curBW = curBW + BWincrement;
        if GainFcprev > GainFc
            d = dprev;
            gainok=1;
        else
            GainFcprev=GainFc;
            dprev =d;
        end
    end
    iterations =iterations +1;
    
end

end




function [d,gainok]=LowPassFIR(BWin,N,Fc,Fs)

gainok =0;
curBW = BWin;
iterations =0;
BWincrement = 25;
GainFcprev =0;
while gainok ==0 && iterations < 10
    
    F6dB2 = Fc+curBW/2;
    
    
    %     h = fdesign.bandpass('n,fc1,fc2', N, F6dB1, F6dB2, Fs);
    %
    %     d = design(h, 'window', 'Window', 'blackmanharris');
    
    d = designfilt('lowpassfir', ...       % Response type
        'FilterOrder',N, ...            % Filter order
        'CutoffFrequency',F6dB2, ...
        'DesignMethod','window', ...         % Design method
        'Window','blackmanharris', ...         % Design method options
        'SampleRate',Fs);               % Sample rate
    
    
    
    [gainok,GainFc]=CheckGainFc(Fc,Fs,d,1e-2);
    
    if ~gainok
        curBW = curBW + BWincrement;
        if GainFcprev > GainFc
            d = dprev;
            gainok=1;
        else
            GainFcprev=GainFc;
            dprev =d;
        end
    end
    iterations =iterations +1;
    
end

end





function [ok,GainFc]=CheckGainFc(Fc,Fs,d,tol)

[H,W] =freqz(d,Fc-10:0.2:Fc+10,Fs);
fc_idx = find (W > Fc,1);
GainFc =abs(H(fc_idx));

ok = 1- GainFc < tol;

end








