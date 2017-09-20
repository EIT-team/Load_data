function [ PhaseAngle ] = ScouseTom_data_PhaseEst( PhaseIn,Protocol,StartInj)
%[PhaseAngle] = ScouseTom_data_PhaseEst( PhaseIn,Protocol,StartInj)
%ScouseTom_data_PhaseEst 
%   Estimate the phase angle on measurement channels by comparing the phase
%   difference to an injection channel. 
%   explanation goes here
%
%   Inputs:
%   PhaseIn - Phase on all channels, output from get_BV
%   Protocol - Complete current injection protocol, used to find which are
%       the injection pairs
%   StartInj - Which line the data starts on
%
%   Outputs: 
%   PhaseAngle - Adjusted phase on all channels (N_prt,N_elec)

%% Read inputs
%get info about data
N_prt=size(Protocol,1);
N_elec=size(PhaseIn,2);
N_inj=size(PhaseIn,1);


if exist('StartInj','var') ==0
    StartInj=1;
    fprintf(2,'StartInj missing from phase est. assuming 1\n');
end



%% Find Phase relative to source elec

%get phase angle for each protocol injection line by subtracting the phase
%from the injection electrode

%make vector of indicies of the protocol. Saves having to calc this each
%time in loop.

%make it longer than necessary to start with
Prt_vecfull=mod(1:N_inj+StartInj+1,N_prt);
Prt_vecfull(Prt_vecfull ==0)=N_prt;

%then take the bit we want
% Prt_vec=Prt_vecfull(StartInj:N_inj+1);
Prt_vec=Prt_vecfull;

PhaseAngleTmp=nan(size(PhaseIn));

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
%% Ensure between [-pi pi]

a = mod(PhaseAngleTmp,2*pi); % [0 2pi)
% shift
j = a > pi;
a(j) = a(j) - 2*pi;
j = a <- pi;
a(j) = a(j) + 2*pi;

PhaseAngle =a;


end

