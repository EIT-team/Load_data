function [ PhaseAngle, PhaseAngleSTD ] = ScouseTom_data_PhaseEst( PhaseIn,Protocol,StartInj)
%ScouseTom_data_PhaseEst Estimate the phase angle through comparison to the
%injection channels
%   Detailed explanation goes here

%get info about data
N_prt=size(Protocol,1);
N_elec=size(PhaseIn,2);
N_inj=size(PhaseIn,1);


if exist('StartInj','var') ==0
    StartInj=1;
    fprintf(2,'StartInj missing from phase est. assuming 1\n');
end

%get phase angle for each protocol injection line by subtracting the phase
%from the injection electrode

%% Find Phase relative to source elec
%make vector of indicies of the protocol. Saves having to calc this each
%time in loop.

%make it longer than necessary to start with
Prt_vecfull=mod(1:N_inj+StartInj+1,N_prt);
Prt_vecfull(Prt_vecfull ==0)=N_prt;

%then take the bit we want
% Prt_vec=Prt_vecfull(StartInj:N_inj+1);
Prt_vec=Prt_vecfull;

PhaseAngleTmp=nan(size(PhaseIn));
STDtmp=PhaseAngleTmp;

% loop through each injection
%and each repeat
for iInj=1:N_inj
    
    %take only the specific injection
    curPA=PhaseIn(iInj,:);
    %and find difference relative to the first line of the injection
    %protocol - this is the source. Also unwrap it to give a more
    %sensible value
    curPA=(curPA - repmat(curPA(Protocol(Prt_vec(iInj),1)),1,N_elec));
    
    %the phase angle is estimated as the mean of the data, within the
    %section of interest
    PhaseAngleTmp(iInj,:)=curPA;
     
end
%% Do some extra stuff...?

%unwrap results now its in the correct format - this isnt needed?
PhaseAngle=unwrap(PhaseAngleTmp);



end

