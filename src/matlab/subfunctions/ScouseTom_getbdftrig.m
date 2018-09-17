function [ StatusChn,TrigPos ] = ScouseTom_getbdftrig( HDR,varargin )
%[StatusChn,TrigPos] = ScouseTom_getbdftrig(HDR)
%SCOUSETOM_GETBDFTRIG Converts bdf triggers into separate channels
%   BioSemi stores status of all 16 channels as 24 bit value, so turn this
%   into a vector of 8 bits. HDR already stores only when things have changed
%
%   Inputs:
%   HDR - Structure from ScouseTom_getHDR
%   Trignum(8) - User specified number of channels to load (optional)
%
%   Outputs:
%   StatusChn  - The state of each channel as boolean at each time point
%   TrigPos     - The time in samples when these changes in digital
%       triggers occured
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

