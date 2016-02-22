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
% it myself BECAUSE IT FUCKING DELETES HALF OF THE CODES FOR NO REASON FOR
% FUCKS SAKE THATS LIEK HOURS OF MY LIFE WASTED WHAT THE
% FUCK!!!!!!
%% Find Marker Codes

%get marker codes using code from BVA import
MRK=readbvconf(HDR.FILE.Path,[HDR.FILE.Name '.vmrk']);
EVENT=parsebvmrk(MRK);

%pull it out of the annoying structure array
latency = [EVENT.latency].';
type = {EVENT.type}.';


 TrigPosOrig=latency; %get the samples at which each of these occur
 [Codes, CodeIndex, Coderef] = unique(type);


% TrigPosOrig=HDR.EVENT.POS; %get the samples at which each of these occur
% [Codes, CodeIndex, Coderef] = unique(HDR.EVENT.Desc);


%% Convert Codes to Decimal

Codes_dec=nan(size(Codes));

%convert each code to decimal number
for iCode=1:length(Codes)
    
    tmp=str2double(Codes{iCode}(2:end));
    
    %if its empty or has some other type e.g. 'boundary'
    if ~isempty(tmp) && ~isnan(tmp)
    Codes_dec(iCode)=str2double(Codes{iCode}(2:end));
    else
        Codes_dec(iCode)=0;
    end
end

%find which ones start with R, and thus need to be shifted 4 bits
R_codes=find(~cellfun(@isempty,(regexpi(Codes,'R'))))';
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

% readbvconf() - read Brain Vision Data Exchange format configuration 
%                file
%
% Usage:
%   >> CONF = readbvconf(pathname, filename);
%
% Inputs:
%   pathname  - path to file
%   filename  - filename
%
% Outputs:
%   CONF      - structure configuration
%
% Author: Andreas Widmann, University of Leipzig, 2007

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2007 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: readbvconf.m 44 2009-11-12 02:00:56Z arnodelorme $
function CONF = readbvconf(pathname, filename)
%this is a function from 


if nargin < 2
    error('Not enough input arguments');
end

% Open and read file
[IN, message] = fopen(fullfile(pathname, filename));
if IN == -1
    [IN, message] = fopen(fullfile(pathname, lower(filename)));
    if IN == -1
        error(message)
    end;
end
raw={};
while ~feof(IN)
    raw = [raw; {fgetl(IN)}];
end
fclose(IN);

% Remove comments and empty lines
raw(strmatch(';', raw)) = [];
raw(cellfun('isempty', raw) == true) = [];

% Find sections
sectionArray = [strmatch('[', raw)' length(raw) + 1];
for iSection = 1:length(sectionArray) - 1

    % Convert section name
    fieldName = lower(char(strread(raw{sectionArray(iSection)}, '[%s', 'delimiter', ']')));
    fieldName(isspace(fieldName) == true) = [];

    % Fill structure with parameter value pairs
    switch fieldName
        case {'commoninfos' 'binaryinfos'}
            for line = sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1
                splitArray = strfind(raw{line}, '=');
                CONF.(fieldName).(lower(raw{line}(1:splitArray(1) - 1))) = raw{line}(splitArray(1) + 1:end);
            end
        case {'channelinfos' 'coordinates' 'markerinfos'}
            for line = sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1
                splitArray = strfind(raw{line}, '=');
                CONF.(fieldName)(str2double(raw{line}(3:splitArray(1) - 1))) = {raw{line}(splitArray(1) + 1:end)};
            end
        case 'comment'
            CONF.(fieldName) = raw(sectionArray(iSection) + 1:sectionArray(iSection + 1) - 1);
    end
end
end

% parsebvmrk() - convert Brain Vision Data Exchange format marker
%                configuration structure to EEGLAB event structure
%
% Usage:
%   >> EVENT = parsebvmrk(MRK);
%
% Inputs:
%   MRK   - marker configuration structure
%
% Outputs:
%   EVENT - EEGLAB event structure
%
% Author: Andreas Widmann, University of Leipzig, 2007

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2007 Andreas Widmann, University of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id: parsebvmrk.m 37 2007-06-26 12:56:17Z andreaswidmann $

function EVENT = parsebvmrk(MRK)

for idx = 1:length(MRK.markerinfos)
    [mrkType mrkDesc EVENT(idx).latency EVENT(idx).duration  EVENT(idx).channel EVENT(idx).bvtime] = ...
        strread(MRK.markerinfos{idx}, '%s%s%f%d%d%d', 'delimiter', ',');

    if strcmpi(mrkType, 'New Segment') || strcmpi(mrkType, 'DC Correction')
        EVENT(idx).type = 'boundary';
    else
        EVENT(idx).type = char(mrkDesc);
    end

    EVENT(idx).code = char(mrkType);
    EVENT(idx).urevent = idx;
end
end


