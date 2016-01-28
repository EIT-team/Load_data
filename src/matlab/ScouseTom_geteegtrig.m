function [ StatusChn,TrigPos ] = ScouseTom_geteegtrig( HDR,varargin )
%SCOUSETOM_TRIG Summary of this function goes here
%   Detailed explanation goes here

%% Check inputs
if ~isempty(varargin)
    
    trignum =varargin{1};
else
    
    trignum=8;
end


% The 8 trigger channels are separated into 1-4 "Stimulus" and 5-8
% "Response". So S3 is 1100 on chn 1-4 and R2 is 0100 on 5-8. BUT these are
% stored in EVENT.TYP as 3 and 2. we want to combine all into a single 8
% bit binary vector. Thus we need to add 16 to the markers for Response
% channel. If a pin is high on both S and R, you get duplicate markers with
% the same time code. so these need to be combined too

% The output in HDR.TYP gives weird results, so I have ignored it and done
% it myself
%% Find Marker Codes

TrigPosOrig=HDR.EVENT.POS; %get the samples at which each of these occur
[Codes, CodeIndex, Coderef] = unique(HDR.EVENT.Desc);


%% Convert Codes to Decimal

Code_dec=nan(size(Codes));

%convert each code to decimal number
for iCode=1:length(Codes)
    
    tmp=str2double(Codes{iCode}(2:end));
    
    if ~isempty(tmp)
    Codes_dec(iCode)=str2double(Codes{iCode}(2:end));
    else
        Codes_dec(iCode)=0;
    end
end

%find which ones start with R, and thus need to be shifted 4 bits
R_codes=find(~cellfun(@isempty,(regexpi(Codes,'R'))));
%shift the relavant ones
Codes_dec(R_codes)=bitshift(Codes_dec(R_codes),4);

%% Find decimal values of each status marker and trim
StatusChnDecFull=Codes_dec(Coderef)';

%remove zero codes as these are not useful to us
rem_idx=find(StatusChnDecFull ==0);
StatusChnDecFull(rem_idx) =[];
TrigPosOrig(rem_idx)=[];

%% Combine duplicate timecodes - this is combining the R and S

%find the unique trigger times, and their index
[TrigPos, ~, idx] = unique(TrigPosOrig);

%Sum up the values for the duplicate array (thanks google). This combines
%the vlaues on the R and S channels to a single value for each time point
StatusChnDec=accumarray(idx,StatusChnDecFull);

%% Make Binary Vector
StatusChn=dec2bin(StatusChnDec,trignum)-'0';%convert into binary vector
StatusChn=fliplr(StatusChn); %sort into LSB


end

