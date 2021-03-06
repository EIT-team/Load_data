function [ Z,Zstd ] = ScouseTom_data_estZ( BV,Elec_inj,ZScaleFactor)
% [ Z,Zstd ] = ScouseTom_data_estZ( BV,Elec_inj,ZScaleFactor)
%SCOUSETOM_DATA_ESTZ Estimates the contact impedance on all electrodes from
%the voltages during injection. Result is average of all injections using
%that electrode. Returns nans if electrodes are unused
%
% Normally called from ProcessBV

%% preallocate

N_freq=size(BV,2);
N_elec=size(Elec_inj,2);

Z=cell(N_freq,1);

for iFreq=1:N_freq
    N_rep=size(BV{iFreq},2);
    Z{iFreq}=nan(N_elec,N_rep);
end
Zstd=Z;

%% Estimate contact impedance

for iFreq=1:N_freq
    
    %get contact impedance values from BV
    for iElec=1:N_elec
        
        %set to nan is there were no injections on this elec
        if isempty(Elec_inj{iElec})
            Z{iFreq}(iElec,:)=nan;
            Zstd{iFreq}(iElec,:)=nan;
        else %average all the voltages on the injection channel
            % Contact impedance is Voltage(V) X 1/Current(A). Zscalefactor
            % corrects for uV and current given in uA
            Z{iFreq}(iElec,:)=ZScaleFactor(iFreq)*nanmean(BV{iFreq}(Elec_inj{iElec},:),1);
            Zstd{iFreq}(iElec,:)=ZScaleFactor(iFreq)*nanstd(BV{iFreq}(Elec_inj{iElec},:),1);
        end
    end
    
end


end

