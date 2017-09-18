function [ trim_demod,FilterOut,Fc] = ScouseTom_data_GetFilterTrim( Vseg,Fs,BWtarget,MaxImpSamples,Fc,plotflag )
% [trim_demod,FilterOut,Fc] = ScouseTom_data_GetFilterTrim( Vseg,Fs,BWtarget,MaxImpSamples,Fc,plotflag )
%
%ScouseTom_GetFilterTrim Gets optimal filter parameters and samples to trim
%   based on SNR tests for different windows.
%

% from the dexterous exploratory hands of jimmy


%% Check inputs

if  exist('Fc','var') == 0 || isempty(Fc)
    %get carrier frequency
    Fc=ScouseTom_data_GetCarrier(Vseg,Fs);
end

if  exist('BWtarget','var') == 0 || isempty(BWtarget)
    BWtarget =100;
end


%get number of samples in segment
Nsamples=length(Vseg);


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

%amount filter ripple must decay by before using data this is based on
%%number of samples in chapter 3 in my thesis. although it could be more
%rigourously chosen
Decay_coef=0.0001;

% Set default parameters for filters. Chosen as trade off between filter
% coefficients and SNR 

%minimum stop band frequency, below this things break
Fstopomin=0.5; 
%minimum frequency at which IIR filters behave - carrier frequencies below
%this need to use FIR regardless
FcominIIR = 15; 
%Stop band attenuation in dB
Astop = 60;
% range of potential transition bands. 
StopBandDiffMax = 350;
StopBandDiffMin = 50;

% minimum frequencyes which have to use FIR bandpass, lower than this its
% best to use separate high and low
FcominFIR =8;


% Frequency below which FIR filters are best separeated into high and low
% pass rather than a single bandpass
FIRLowPassCutoff = 20;
% equivalent for IIR
IIRLowPassCutoff = 16;

%flags for using low pass only
UseLowPassIIR=0;
UseLowPassFIR=0;


%the filter takes some time in seconds to decay to ~zero for IIR butterworth filters.
%This is independent of sampling rate, so we need to have a different max
%number of samples for different sampling rates. This is also the max order
%for FIR. This was chosen based on chap3 of my thesis, not super robustly

Mindecay_seconds=0.2;
MaxFirOrder = floor(Mindecay_seconds*Fs);


%% Find best IIR filter
% We want to use IIR when we can but it has problems: it can go unstable,
% and take too long to decay. we aim for having ~50% of the injection time
% unafected by filtering artefacts. So we need to check the IIR filter
% impulse response

%set flags 
IIRfailed =0;
dprev=[];

try
    %start with target bandwidth, this will increase if not
    BWIIR = BWtarget;
    IIRfound =0;
    curStopBandDiff = StopBandDiffMax;
    
    while ~IIRfound
        
        % find the IIR filter for this FC, either bandpass or
        % lowpass depending on the FC
        
        % find the best filter for this given transition band
        if Fc > IIRLowPassCutoff
            [dcur,st,gainok,implen]=BandPassIIR(BWIIR,curStopBandDiff,Astop,Fc,FcominIIR,Fstopomin,Fs,Decay_coef);
        else
            UseLowPassIIR=1;
            [dcur,st,gainok,implen]=LowPassIIR(BWIIR,curStopBandDiff,Astop,Fc,Fs,Decay_coef);
        end
        
        %if the filter is both stable and below the maximum measurement
        %time, then we can try reducing the transition band for better SNR
        if st && implen < MaxImpSamples
            
            curStopBandDiff = curStopBandDiff - 50;
            dprev = dcur;
            
            % if  we have reached the minimum transition time, or we have
            % broken the gain, then stop
            if curStopBandDiff < StopBandDiffMin || ~gainok
                IIRfound =1;
            end
            
        else
            % if the impulse response is greater than maxsamples,  then
            % stop and use the last one
            IIRfound =1;
            
        end
        
        
    end
    
    % the best possible IIR from the previous loop
    dIIR = dprev;
    
    %if it broke then set how many samples we need to trim
    if isempty(dIIR)
        IIRfailed =1;
        Samples_needed = inf;
        IIRgainok =0;
    else
        
        Samples_needed = impzlength(dIIR,Decay_coef);
        
    end
    
    % plot impulse response if required, indicating max samples and minimum
    % samples needed for filter to decay
    if plotflag ==1 && ~IIRfailed
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
% prevent ridiculous orders of filters
N = min([MaxImpSamples MaxFirOrder]);

% create filter with narrowest bandwidth possible, either low or bandpass
% given FC. Blackman harris window chosen as it gives the best trade off
% between stopband ripple and rolloff
if Fc > FIRLowPassCutoff
    [dFIR,FIRgainok]=BandPassFIR(BWFIR,N,Fc,FcominFIR,Fs);
else
    [dFIR,FIRgainok]=LowPassFIR(BWFIR,N,Fc,Fs);
    UseLowPassFIR=1;
end
%% choose filter

% sub 1Hz carriers must use FIR anyway
if Fc < 1
    ForceFIR =1;
else
    ForceFIR =0;
end

%from TestFilterSNR - FIR outperforms IIR until around 1200 samples (on biosemi Fs and BW50), then
%they are within 1% of each other. IIR is MUCH faster than high order FIR
%so use this to speed up the process

if (MaxImpSamples < Samples_needed) || IIRfailed || ForceFIR
    %if we do not have enough samples to allow for the filter to decay
    %sufficiently, then use the slower FIR filter. 
    
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
%% plot impulse and frequency response
if plotflag ==1
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
% Creates the bandpass IIR filters. Checks the gain at the carrier
% frequency is unaffected. If this is the case it increases the bandwidth.
% Checks for stability and calculates the time for filter time response to
% decay to zero.

Apass  = 0.5;

gainok =0;
curBW = BWin;
iterations =0;
BWincrement = 25;
GainFcprev =0;

maxiterations =10;

while gainok ==0 && iterations < maxiterations
    
    % lowest cuttoff
    Fpass1bp = Fc-curBW/2;
    
    %stopnand attenuation
    Astop1 = Astop;
    Astop2= Astop;
    
    if (Fpass1bp <= Fcomin)
        Fpass1bp = max([Fc/2 Fcomin]);
        Astop1 = 15;
    end
    
    % lowest stop band freq
    Fstop1bp = Fpass1bp -StopBandDiff;
    Fstop1bp = max([Fstop1bp Fstopomin]);
    
    %highest freq cut off
    Fpass2bp = Fc+curBW/2;
    Fstop2bp = Fpass2bp + StopBandDiff;
    
    %build filter
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
    
    
    st = isstable(d); % check the stability of this filter
    implen=impzlength(d,Decay_coef); %how long does it take to decay?
    
    %check the gain
    [gainok,GainFc]=CheckGainFc(Fc,Fs,d,1e-2);
    
    %if the gain is ok then change the bandwidth
    if ~gainok
        curBW = curBW + BWincrement;
        %if we have made things worse then stop
        if GainFcprev > GainFc
            d = dprev;
            gainok=0;
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
% Creates the lowpass IIR filters. Checks the gain at the carrier
% frequency is unaffected. If this is the case it increases the bandwidth.
% Checks for stability and calculates the time for filter time response to
% decay to zero.


Apass  = 0.5;

gainok =0;
curBW = BWin;
iterations =0;
BWincrement = 25;
GainFcprev =0;
while gainok ==0 && iterations < 10
    
    Fpass = Fc+curBW/2;
    Fstop= Fpass + StopBandDiff;
    
    
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
    
    %if the gain is not ok, then increase the cutoff frequency and try
    %again
    if ~gainok
        curBW = curBW + BWincrement;
        % if we have made things worse, then stop
        if GainFcprev > GainFc
            d = dprev;
            gainok=1;
        else
            % else move on to the next iteration
            GainFcprev=GainFc;
            dprev =d;
        end
    end
    iterations =iterations +1;
    
end

end

function [d]=HighPassIIR(Fc,Fs)
% makes high pass IIR. This is always the same, so no iterations for this

Fpass = max([(Fc -10) 1.5]);
Fstop = max([(Fpass -1.5) 1]);
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
% Creates a bandpass FIR filter, checking the gain at the carrier is
% unaffected. Increases the bandwidth until this is the case
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
    
    %build filter
    d = designfilt('bandpassfir', ...       % Response type
        'FilterOrder',N, ...            % Filter order
        'CutoffFrequency1',F6dB1, ...    % Frequency constraints
        'CutoffFrequency2',F6dB2, ...
        'DesignMethod','window', ...         % Design method
        'Window','blackmanharris', ...         % Design method options
        'SampleRate',Fs);               % Sample rate
    
    
    % check gain
    [gainok,GainFc]=CheckGainFc(Fc,Fs,d,1e-2);
    
    %if the gain is not ok, then increase the cutoff frequency and try
    %again
    if ~gainok
        curBW = curBW + BWincrement;
        % if we have made things worse, then stop
        if GainFcprev > GainFc
            d = dprev;
            gainok=1;
        else
            % else move on to the next iteration
            GainFcprev=GainFc;
            dprev =d;
        end
    end
    iterations =iterations +1;
    
end

end




function [d,gainok]=LowPassFIR(BWin,N,Fc,Fs)
%creates a low pass FIR filter

gainok =0;
curBW = BWin;
iterations =0;
BWincrement = 25;
GainFcprev =0;
while gainok ==0 && iterations < 10
    
    % set cutoff freq
    F6dB2 = Fc+curBW/2;
    
    % build filter
    d = designfilt('lowpassfir', ...       % Response type
        'FilterOrder',N, ...            % Filter order
        'CutoffFrequency',F6dB2, ...
        'DesignMethod','window', ...         % Design method
        'Window','blackmanharris', ...         % Design method options
        'SampleRate',Fs);               % Sample rate
    
    
    % check if gain is ok
    [gainok,GainFc]=CheckGainFc(Fc,Fs,d,1e-2);
    
    %if the gain is not ok, then increase the cutoff frequency and try
    %again
    if ~gainok
        curBW = curBW + BWincrement;
        % if we have made things worse, then stop
        if GainFcprev > GainFc
            d = dprev;
            gainok=1;
        else
            % else move on to the next iteration
            GainFcprev=GainFc;
            dprev =d;
        end
    end
    iterations =iterations +1;
    
end

end

function [ok,GainFc]=CheckGainFc(Fc,Fs,d,tol)
% checks if the gain at the carrier frequecy is altered by the filtering,
% within a given tolerance

[H,W] =freqz(d,Fc-10:0.2:Fc+10,Fs);
fc_idx = find (W > Fc,1);
GainFc =abs(H(fc_idx));

ok = 1- GainFc < tol;

end

