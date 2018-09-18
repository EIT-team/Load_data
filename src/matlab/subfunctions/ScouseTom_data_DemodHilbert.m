function [ Vdata_demod,Pdata_demod,data ] = ScouseTom_data_DemodHilbert( data,Filt)
% [Vdata_demod,Pdata_demod,data] = ScouseTom_data_DemodHilbert( data,Filt)
%demod_hilbert - filters and demodulates data using hilbert transform
%method.
%   Inputs:
%   Data- signal to be demodulated
%   Filt - Filter object


%% Check data
if any(isnan(data))
    Vdata_demod = nan;
    Pdata_demod=nan;
    fprintf(2,'Nans in data!\n');
    return
end

%% Filter data

% remove DC offset
datamean = mean(data);
data = bsxfun(@minus,data,datamean);

% multiple filters can be specified in cell array
if ~iscell(Filt)
    %filter data using coefs
    data = filtfilt(Filt,data);
else
    %apply in filter in turn
    for iFilter = 1 : length(Filt)
        data = filtfilt(Filt{iFilter},data);
    end
    
end
%% Demodulate


% remove DC offset
datamean = mean(data);
data = bsxfun(@minus,data,datamean);


%get envelope of signal using hilbert transform
data = (hilbert(data));

Vdata_demod=abs(data); %amplitude is abs of hilbert
Pdata_demod=angle(data); % phase is the imaginary part



end