function [ StatusChn,TrigPos ] = ScouseTom_getbdftrig( HDR,varargin )
%SCOUSETOM_GETBDFTRIG Converts bdf triggers into separate channels
%  User can specify the number of trigger channels to load, but its
%  basically always 8
%% Check inputs
if ~isempty(varargin)
    
    trignum =varargin{1};
else
    
    trignum=8;
end

%% Convert 24bit value into binary vector
StatusChn=dec2bin(HDR.BDF.Trigger.TYP)-'0';%convert into binary vector
StatusChn=StatusChn(:,end-(trignum-1):end); % take only last 8 bits
StatusChn=fliplr(StatusChn); %sort into LSB
TrigPos=HDR.BDF.Trigger.POS; %get the samples at which each of these occur


end

