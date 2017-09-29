function [Vmag,Vphase,Vmag_std,Vphase_std]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod,InjectionWindowsFull)
% [Vmag,Vphase,Vmag_std,Vphase_std]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,Trim_demod,InjectionWindowsFull)
%
%ScouseTom_data_getBV Takes the average of the demodulates voltages within
%each injection. This gives the conventional "boundary voltage" used for
%EIT reconstructions. 

Nsample=size(Vdata_demod,1);
Nelec=size(Vdata_demod,2);

% Find switches within the data given
CompleteInj=all(InjectionWindowsFull < Nsample,2);

% Extract the relevant injection windows
InjWind=InjectionWindowsFull(CompleteInj,:);
% trim off edge artefacts from filtering
InjWind(:,1)=InjWind(:,1)+Trim_demod;
InjWind(:,2)=InjWind(:,2)-Trim_demod;

nInj=size(InjWind,1);

%% preallocate
Vmag = nan(nInj,Nelec);
Vphase=nan(size(Vmag));
Vmag_std=nan(size(Vmag));
Vphase_std=nan(size(Vmag));

%% Loop through each injection window
for iInj = 1:nInj
    
    %get the current window
    curWind=InjWind(iInj,:);
    
    %take the values of interest
    curV= Vdata_demod(curWind(1):curWind(2),:);
    curP=unwrap(Pdata_demod(curWind(1):curWind(2),:));
    
    %magnitude
    Vmag(iInj,:)=nanmean(curV);
    Vmag_std(iInj,:)=nanstd(curV);
    %phase
    Vphase(iInj,:)=nanmean(curP);
    Vphase_std(iInj,:)=nanstd(curP); 
end

end