function [ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( data,B,A)
%demod_hilbert - filters and demodulates data using hilbert transform
%method. Slight modification of G-Dragons code get_BV2
%   Inputs:
%   Data- signal to be demodulated
% B A - Filter coefficients

if any(isnan(data))
    Vdata_demod = nan;
    Pdata_demod=nan;
    
    fprintf(2,'Nans in data!');
    return
end


if ~iscell(A)
    
    %filter data using coefs
    data = filtfilt(B,A,data);
    
else
    for iFilter = 1 : length(A)
        data = filtfilt(B{iFilter},A{iFilter},data);
    end
    
end

%get envelope of signal using hilbert transform
data = (hilbert(data));

Vdata_demod=abs(data); %amplitude is abs of hilbert
Pdata_demod=angle(data); % phase is the imaginary part



end

