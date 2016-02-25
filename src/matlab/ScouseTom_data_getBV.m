function     [Vmag,Phase]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod,InjectionWindowsFull)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nsample=size(Vdata_demod,1);

CompleteInj=all(InjectionWindowsFull > Nsample);

InjWind=InjectionWindowsFull(CompleteInj,:);

InjWind(:,1)=InjWind(:,1)+Trim_demod;
InjWind(:,2)=InjWind(:,2)-Trim_demod;

nInj=size(InjWind,1);

%% preallocate


Vmag = nan(

%%
for iInj = 1:nInj
    
    curWind=InjWind(iInj,:);
    
    curV= Vdata_demod(curWind(1):curWind(2));
    curP=Pdata_demod(curWind(1):curWind(2));
    




end

