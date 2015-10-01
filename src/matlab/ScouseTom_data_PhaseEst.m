function [ PhaseAngle, PhaseAngleSTD ] = ScouseTom_data_PhaseEst( Pdemod,trim_demod,Protocol,varargin )
%ScouseTom_data_PhaseEst Estimate the phase angle through comparison to the
%injection channels
%   Detailed explanation goes here

%get info about data
N_prt=size(Pdemod,1);
N_elec=size(Pdemod,3);
N_rep=size(Pdemod,4);

%get phase angle for each protocol injection line by subtracting the phase
%from the injection electrode - ignoring trim section

%%

PhaseAngleTmp=nan(N_prt,N_elec,N_rep);
STDtmp=PhaseAngleTmp;

% loop through each protcol line
for iPrt=1:N_prt
%and each repeat
  for iRep=1:N_rep  
    
      %take only the specific injection
      curPA=(squeeze(Pdemod(iPrt,:,:,iRep)));
      %and find difference relative to the first line of the injection
      %protocol - this is the source. Also unwrap it to give a more
      %sensible value 
      curPA=unwrap(curPA - repmat(curPA(:,Protocol(iPrt,1)),1,N_elec));
      
      %the phase angle is estimated as the mean of the data, within the
      %section of interest
    PhaseAngleTmp(iPrt,:,iRep)=(nanmean(curPA(trim_demod:end-trim_demod,:),1));
    %also store the standard dev, as an indication of the quality of the
    %estimate
    PhaseAngleSTDTmp(iPrt,:,iRep)=nanstd(curPA(trim_demod:end-trim_demod,:),1);

  end
end
%%
if (isempty(varargin) == 1 && N_prt > 1)%if no extra input then go ahead and reshape into the nromal format of Prot x Rep
    
    %put each protcol line one after each other
    for iPrt =1:N_prt
        idx=((iPrt-1)*N_elec) +1;
        
        if N_rep ==1
            PhaseAngle(idx:idx+N_elec-1,1)=PhaseAngleTmp(iPrt,:);
            PhaseAngleSTD(idx:idx+N_elec-1,1)=PhaseAngleSTDTmp(iPrt,:);
            
        else

        PhaseAngle(idx:idx+N_elec-1,:)=squeeze(PhaseAngleTmp(iPrt,:,:));
        PhaseAngleSTD(idx:idx+N_elec-1,:)=squeeze(PhaseAngleSTDTmp(iPrt,:,:));
        end
        
    end
    
    %unwrap results now its in the correct format
    PhaseAngle=unwrap(PhaseAngle,[],2);
    
    
    
else %i dont know why you would need this
    PhaseAngle=PhaseAngleTmp;
    PhaseAngleSTD=PhaseAngleSTDTmp;
end


end

