function [ Vdata_demod,Pdata_demod,data ] = ScouseTom_data_DemodHilbert( data,Filt)
%demod_hilbert - filters and demodulates data using hilbert transform
%method. Slight modification of G-Dragons code get_BV2
%   Inputs:
%   Data- signal to be demodulated
%   Filt - Filter object

if any(isnan(data))
    Vdata_demod = nan;
    Pdata_demod=nan;
    
    fprintf(2,'Nans in data!\n');
    return
end

%%

if ~iscell(Filt)
    
    %filter data using coefs
    data = filtfilt(Filt,data);
    
else
    for iFilter = 1 : length(Filt)
        data = filtfilt(Filt{iFilter},data);
    end
    
end

%get envelope of signal using hilbert transform
data = (hilbert(data));

Vdata_demod=abs(data); %amplitude is abs of hilbert
Pdata_demod=angle(data); % phase is the imaginary part



end
