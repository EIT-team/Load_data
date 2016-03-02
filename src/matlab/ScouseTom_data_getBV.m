function [Vmag,Vphase,Vmag_std,Vphase_std]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod,InjectionWindowsFull)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nsample=size(Vdata_demod,1);

CompleteInj=all(InjectionWindowsFull < Nsample,2);

InjWind=InjectionWindowsFull(CompleteInj,:);

InjWind(:,1)=InjWind(:,1)+Trim_demod;
InjWind(:,2)=InjWind(:,2)-Trim_demod;

nInj=size(InjWind,1);

%% preallocate


Vmag = nan(nInj,1);
Vphase=nan(size(Vmag));
Vmag_std=nan(size(Vmag));
Vphase_std=nan(size(Vmag));

%% Loop through each injection window
for iInj = 1:nInj
    
    %get the current window
    curWind=InjWind(iInj,:);
    
    %take the values of interest
    curV= Vdata_demod(curWind(1):curWind(2));
    curP=Pdata_demod(curWind(1):curWind(2));
    
    %voltage
    Vmag(iInj)=nanmean(curV);
    Vmag_std(iInj)=nanstd(curV);
    
    Vphase(iInj)=nanmean(curP);
    Vphase_std(iInj)=nanstd(curP);
    
    
end

end