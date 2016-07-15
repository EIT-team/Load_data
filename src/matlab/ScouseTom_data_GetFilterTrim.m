function [ trim_demod, B,A,Fc ] = ScouseTom_data_GetFilterTrim( Vseg,Fs,BW,plotflag )
%ScouseTom_GetFilterTrim Gets optimal filter parameters and samples to trim
%   based on SNR tests for different windows
% from the dexterous exploratory hands of jimmy

%get carrier frequency
Fc=ScouseTom_data_GetCarrier(Vseg,Fs);

%get number of samples in segment
Nsamples=length(Vseg);

%trim max is 25% of signal - to give 50% for the average
trim_max=ceil(0.25*Nsamples);
%round up to nearest 10 samples - this is so that data with window size 1
%sample different (as can happen) produces the same filter, and thus the
%same result.
rndnum=10;
trim_max=round(trim_max/rndnum)*rndnum;

%amount filter ripple must decay by before using data%this is based on
%number of samples in chapter 3 in y thesis. although it could be more
%rigourously chosen
Decay_coef=0.001;

%starting order of butterworth filter
Nord=3;

%flag for whether we need to add a separate high pass filter, if things are
%exploding due to narrow bandwidth
AddSecondFilter=0;

%the filter takes some time in seconds to decay to ~zero for IIR butterworth filters.
%This is independent of sampling rate, so we need to have a different max
%number of samples for different sampling rates

decay_seconds=0.08;
threshold_samples = floor(decay_seconds*Fs);

% create time vector for investigating the impulse response of the filter,
% preventing it from being huge for long data sets
if Nsamples > 5*Fs
    tdec=ceil(5*Fs);
else
    tdec=5*fix(Nsamples);
end


%% Estimate decay time for IIR filter
% if IIR decay time this is too long then we have to use much slow FIR
% filters to avoid artefacts

IIRfailed =0;

try
    
    if (Fc - BW/2) > 0
        [B,A] = butter(Nord,(Fc+[-BW/2,BW/2])./(Fs/2));
    else
        fprintf(2,'WARNING! CHOSEN BANDWIDTH TOO LARGE FOR CARRIER FREQUENCY! USING TWO FILTERS INSTEAD\n');
        % maximum bandwidth - with at least 1Hz clearance
        maxBWhalf = floor(Fc)-1;
        
        % try third order low pass filter
        [B,A] = butter(Nord,(Fc+[BW/2])./(Fs/2));
        
        % if max occurs in the second half of the impulse response, the filter
        % is unstable
        [h,t]=impz(B,A,tdec);
        [mm, midx]=max(abs(h));
        
        if midx > length(h)/2
            Nord=2;
            [B,A] = butter(Nord,(Fc+[-maxBWhalf,maxBWhalf])./(Fs/2));
        end
        
%         AddSecondFilter=1;
        
    end
    
    %% check decay of the low pass filter
    %get filter response and estimate when to chop data
    [H, T]= impz(B,A,tdec);
    [maxh, ih]=max(abs(H));
    
    %linear fit of the exponetial decay - when it reaches certain percentage of max
    Htofit=log(abs(H(ih:end)./maxh));
    
    minh=Htofit(end);
    
    iminh=find( Htofit < 0.95*minh,1);
    
    % figure;plot(T(ih:end),Htofit)
    P=polyfit(T(ih:iminh+ih-1),Htofit(1:iminh),1);
    %sample where decay reaches the desired value
    Decay_samples=ceil(log(Decay_coef)/P(1));
    
    %total amount to trim is decay time plus time to max
    Samples_needed=Decay_samples+ih;
    
    if plotflag ==1;
        figure;
        hold on
        plot(T,(H))
        %line([0 length(T)],[maxh*Decay_coef maxh*Decay_coef],'color','r');
        line([Samples_needed Samples_needed],[min(H) max(H)],'color','r');
        line([trim_max trim_max],[min(H) max(H)],'color',[0 0.5 0]);
        hold off
        title('impulse response of IIR filter')
        xlim([0 (ceil(Samples_needed/1000))*1000])
        % set(gca,'Yscale','log');
        legend('Filter response','Req. Trim Samples','Max Trim Samples')
        drawnow
    end
    %% add high pass filter coeffs on it
    
catch
    IIRfailed =1 ;
end


%% choose filter

%from TestFilterSNR - FIR outperforms IIR until around 1200 samples (on biosemi Fs and BW50), then
%they are within 1% of each other. IIR is MUCH faster than high order FIR
%so use this to speed up the process

if (trim_max <Samples_needed) || IIRfailed
    %if we do not have enough samples to allow for the filter to decay
    %sufficiently, then use the slower FIR filter. Blackman harris window
    %chosen as it gives the best trade off between stopband ripple and
    %rolloff
    
    filter_size = min([trim_max threshold_samples]);
    
    if (Fc - BW/2) > 0
        [B,A]=fir1(filter_size,(Fc+[-BW/2,BW/2])./(Fs/2),'bandpass',window(@blackmanharris,filter_size+1));
    else
        % maximum bandwidth with 3Hz clearance for slower rolloff
        maxBWhalf = floor(Fc)-3;
        [B,A]=fir1(filter_size,(Fc+[BW/2])./(Fs/2),'lowpass',window(@blackmanharris,filter_size+1));
        fprintf(2,'WARNING! CHOSEN BANDWIDTH TOO LARGE FOR CARRIER FREQUENCY! USING TWO FILTERS\n');
    end
    
    trim_demod=filter_size;
    
    disp('FIR with Blackman-Harris Window used');
else
    %if we have more samples than we need, still use the max as we want the
    %filter to decay as much as possible (I think)
    trim_demod=trim_max;
    
    disp('3rd Order Butterworth Filter Used');
    
end
%%
if plotflag ==1;
    figure;
    impz(B,A);
    drawnow
    figure;
    freqz(B,A);
    drawnow
end

%%

if AddSecondFilter
    B={B};
    A={A};
    [B{2},A{2}] = butter(1,(Fc+[-maxBWhalf])./(Fs/2),'high');
end






end

