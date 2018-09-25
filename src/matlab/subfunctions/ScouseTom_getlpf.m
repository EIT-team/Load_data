function [Filtout,TrimSamples] = ScouseTom_getlpf(n,Fc,Fs)
% [Filtout] = ScouseTom_getlpf(n,Fc,Fs)
% makes butterworth lowpass filter object of specified order
%
% inputs:
% n - filter order
% Fc - Cut of frequencies in Hz i.e. [1000 2000]
% Fs - Sample rate
%
% Output
% Filtout - object so you can use it like filtfilt(Filtout,Data)
% TrimSamples - number of samples contaminated by fltering artefacts you
% need to remove when using this filter


% how much impulse response has to decay by until we trust the data. This
% is set very high for parallel stuff
Decay_coef=1e-9;

Filtout = designfilt('lowpassiir', 'FilterOrder', n, ...
                     'HalfPowerFrequency', Fc, 'SampleRate', Fs*1, ...
                     'DesignMethod', 'butter');

TrimSamples=impzlength(Filtout,Decay_coef);

end


