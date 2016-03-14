function [ Inj_pairs, badnessflag,RMS_CHN ] = ScouseTom_data_EstInjPair( Vin,Threshold )
%ScouseTom_EstimateInjPair Simply estimates which channels are injecting by
%finding the two biggest channels using RMS
%   Detailed explanation goes here

%% parse inputs

if exist('Threshold','var')==0
    Threshold=0.5;
end

N_elec=size(Vin,2);
RMS_CHN=nan(N_elec,1);
badnessflag=0;

%% Find RMS

%get rms for each channel
for iElec=1:N_elec
    tmp=Vin(:,iElec); %only take first bit to save processing time on big injects
    RMS_CHN(iElec)=sqrt(mean(detrend(tmp).^2));
end

%% Find Biggest Two Channels
%sort into order
[RMS_CHNS ,I]=sort(RMS_CHN,1,'descend');


%if the biggest two arent more than 50% than the next one then warn
if RMS_CHNS(3) > Threshold*RMS_CHNS(2)
    fprintf(2,'Voltages on other channels greater than %.2f the voltage on injection pairs\n',Threshold);
    badnessflag=1;
end

Inj_pairs=sort(I(1:2));



end

