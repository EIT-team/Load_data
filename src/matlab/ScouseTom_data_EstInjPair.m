function [ Inj_pairs, badnessflag ] = ScouseTom_data_EstInjPair( Vin )
%ScouseTom_EstimateInjPair Simply estimates which channels are injecting by
%finding the two biggest channels using RMS
%   Detailed explanation goes here

%%

N_elec=size(Vin,2);
RMS_CHN=nan(N_elec,1);
badnessflag=0;

maxlen=1000;

%%

if size(Vin,1) > maxlen
    idxmax=maxlen;
else
    idxmax=size(Vin,1);
end

%get rms for each channel
for iElec=1:N_elec
    tmp=Vin(1:idxmax,iElec); %only take first bit to save processing time on big injects
    RMS_CHN(iElec)=sqrt(mean(detrend(tmp).^2));
end

%sort into order
[RMS_CHN ,I]=sort(RMS_CHN,1,'descend');


%if the biggest two arent more than 50% than the next one then warn
if RMS_CHN(3) > 0.5*RMS_CHN(2)
    warning('Voltages on other channels greater than half the voltage on injection pairs');
    badnessflag=1;
end

Inj_pairs=sort(I(1:2));



end

