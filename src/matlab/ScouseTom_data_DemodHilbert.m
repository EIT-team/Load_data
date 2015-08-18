function [ data_demod ] = ScouseTom_data_DemodHilbert( data,B,A)
%demod_hilbert - filters and demodulates data using hilbert transform
%method. Slight modification of G-Dragons code get_BV2
%   Inputs:
%   Data- signal to be demodulated
% B A - Filter coefficients

% ADD PHASE ESTIMATE INTO THIS BIT!!

if any(isnan(data))
    data_demod = nan;
    
    %     warning('Nans in data bro!');
else
    
    
    %filter data using coefs
    data_demod = filtfilt(B,A,data);
    %get envelope of signal using hilbert transform
    data_demod = abs(hilbert(data_demod));
    
    
    
end

end

