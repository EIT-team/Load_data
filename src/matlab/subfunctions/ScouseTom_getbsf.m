function [Filtout,TrimSamples] = ScouseTom_getbsf(n,Fc,Fs)
% [Filtout] = ScouseTom_getbsf(n,Fc,Fs)
% makes butterworth bandstop filter object of specified order
%
% inputs:
% n - filter order
% Fc - array of cut of frequencies in Hz i.e. [1000 2000]
% Fs - Sample rate
%
% Output
% Filtout - object so you can use it like filtfilt(Filtout,Data)
% TrimSamples - number of samples contaminated by fltering artefacts you
% need to remove when using this filter


% how much impulse response has to decay by until we trust the data. This
% is set very high for parallel stuff
Decay_coef=1e-9;

% using this object allows for higher orders easily
Filtout =designfilt('bandstopiir', 'FilterOrder', n, ...
    'HalfpowerFrequency1',Fc(1) ,...
    'HalfpowerFrequency2',Fc(2),...
    'DesignMethod','butter',...
    'SampleRate', Fs);

TrimSamples=impzlength(Filtout,Decay_coef);



end

