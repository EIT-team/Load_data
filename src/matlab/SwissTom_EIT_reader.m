% Copyright (c) 2013 Swisstom AG, Landquart, Switzerland.
% This file is part of EIT-Pioneer-Set. 
% EIT-Pioneer-Set is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version. 
% EIT-Pioneer-Set is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with EIT-Pioneer-Set.  If not, see <http://www.gnu.org/licenses/>.

function out = ST_EIT_reader(fname)
%ST_EIT_reader Reads a Swisstom .eit file.
%  OUT=ST_EIT_reader(FNAME) reads the file named FNAME
% 
%  INPUTS:
%    FNAME  - filename (string)
%
%  OUTPUTS:
%    OUT    - a struct containg all the data in file FNAME.

%             For Format Version = 3, the fields are:
%        formatVersion: 1
%     initialTimestamp: 1.3301e+12
%       finalTimestamp: 1.3301e+12
%       numberOfFrames: 417
%          dataFileName: [1x100 char]
%     dataFileCondition: [1x300 char]
%      dataFileComments: [1x600 char]
%                  tic: [1x1 struct]              
%                                imageRate: 10
%                         injectionCurrent: 5
%                              currentFreq: 100
%                       switchSettlingTime: 20
%                         injectionPattern: 4
%                                 num_elec: 32
%                  sbc: [1x1 struct]
%                                nco_freq: 5
%                                dac_gain: 100
%                                  dac_fs: 1
%                                    idta: 'INNNNNNNNNNNNNNNNNNNNNNNNNNGNNNN'
%                                    mdta: '1NNNNNNNNNNNNNNNNNNNNNNNNNN2NNNN'
%                                   tickl: 0
%                                   tmckl: 0
%                                    nppw: 10
%                                  pga0_g: 1
%                                  pga1_m: 0
%                                  pga1_g: 0
%                                  pga1_o: 0
%                                    nign: 2
%                                    tign: 0
%                                    nsmp: 320
%                                    tsmp: 1
%                                      tr: 0
%                                     ni0: 4
%                                     nq0: 4
%                                     ni1: 4
%                                     nq1: 4
%                                    nraw: 0
%                                   raw_m: 0
%                                  endian: 0
%               frames: [1x417 struct]
%                        absoluteTimestamp: 1.3414e+12
%                                     size: 8488
%                      StartIdentification: 0x53776973   (typecast to hex)
%                        EndIdentification: 0x73746f6d
%                                     type: 0
%                                     size: 8488
%                                timestamp: 2355219
%                                   status: 0
%                                ErrorCode: 0
%                              ecgDataSize: 0
%                          stretchDataSize: 0
%                             bendDataSize: 0
%                         positionDataSize: 0
%              ScanningPatternCodeDataSize: 0  (instead of ScanningPatternCode)
%                 voltageInjectionDataSize: 256
%                               iqDataSize: 8192
%                         nExternalSignals: 0
%                               ecgPayload: []
%                           stretchPayload: []
%                              bendPayload: []
%                          positionPayload: []
%               ScanningPatternCodePayload: []
%                  voltageInjectionPayload: [64x1 double]
%                                iqPayload: [2048x1 double]
%                                      ext: []
%
%  in frames-payload:
%  in External-Signals: Adjusted size: size of bytes includes 
%      code+size+payload, so divide it by 4 minus 2 to get payload-size

% Author : beat mueller <bmu@swisstom.com>
% Project: pioneer set
% Copyright 2013, Swisstom AG
%
%
% $Id: ST_EIT_reader.m 285 2013-07-01 08:46:03Z bmu $

%%% Main Body
[fid msg]= fopen(fname,'r','ieee-be','UTF-8');
if fid==-1, error(msg); end
%whatever happens, close the file
try
   format_version = fread(fid,1,'int32','ieee-be');
   switch(format_version)
       case 3
         out = read_file(fid,3);
      otherwise
         out = readRawData(fid,format_version);
   end
catch err
   fclose(fid);
   rethrow(err);
end
rest = fread(fid);
if ~isempty(rest)
   warning('ST_EIT_reader:incompleteRead','There was more data in the file than was read');
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = readRawData(fid,first_val)
% reads rawIQ data files. The file contains nothing but integers in sets of
% 2048 (1024 IQ measurements)
% the first value has already been read in, we have to incorporate it
% we want to output a normal eit file. We'll cheat and load a template with
% the standard fields
out = [];                                 % just so Matlab doesn't complain
out.formatVersion = 0;                    % mark this file as funny
frame_tmpl = out.frames(1);
c = 1;
first = 1;
while ~feof(fid)
   frame = frame_tmpl;
   if first
      frame.iqPayload = [first_val; fread(fid,2047,'int32','ieee-be')];
      first = 0;
   else
      frame.iqPayload = fread(fid,2048,'int32','ieee-be');
      if isempty(frame.iqPayload);break; end;
   end
   % set dummy voltage data (like TIC)
   frame.voltageInjectionPayload = frame.iqPayload(1:64);  
   out.frames(c) = frame;
   c = c + 1;
end
if length(out.frames(end).iqPayload) ~= 2048
   out.frames(end) = [];
end
out.numberOfFrames = length(out.frames);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = read_file(fid, ver)
% reads EIT file format v1
out.type = 'eit';
switch ver
   case 3
      var = filevarsV3;
   otherwise
      error('ST_EIT_reader:UnkownVersion','Version number not supported');
end
out.formatVersion = ver;

%%% Read general data
for i = 1:size(var,1)
   if strcmp(var{i,2},'utf8')
      eval(['out.' var{i,1} ' = readUTF(fid);']);
   elseif strcmp(var{i,2}(1:5),'chars')
      if length(var{i,2}) > 5
         ln = var{i,2}(6:end);
         if strcmp(ln,'NE')
            ln = num2str(out.tic.num_elec); % must know in advance
         end
         eval(['out.' var{i,1} ' = readChars(fid,' ln ');']);
      else
         eval(['out.' var{i,1} ' = readChars(fid);']);
      end
   else
      eval(['out.' var{i,1} ' = fread(fid, 1,''' var{i,2} ''',''ieee-be'');']);
   end
end

%%% Read frames
for i = 1:out.numberOfFrames
   out.frames(i) = read_frame(fid,ver);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function frame = read_frame(fid,ver)
frame.type = 'frame';
frame.absoluteTimestamp = fread(fid,1,'int64');

frame.size = fread(fid,1,'int32');

%%% Read the header info
switch ver
   case 3
      var = framevarsV3;
   otherwise
      error('ST_EIT_reader:UnkownVersion','Version number not supported');
end
for i = 1:size(var,1)
   eval(['frame.' var{i,1} ' = fread(fid, 1,''' var{i,2} ''',''' var{i,3} ''');']);
end
switch ver
    case 3
        payloads = {'ecgPayload','stretchPayload','bendPayload' ... 
            'positionPayload', 'scanningPatternCodePayload' ...
            'voltageInjectionPayload' 'iqPayload' };
   otherwise
      error('ST_EIT_reader:UnkownVersion','Version number not supported');
end
%%% Read the payloads
for i = 1:numel(payloads)
   % divide DataSize by four (the size of int);
   frame.(payloads{i}) = ...
      fread(fid,frame.(strrep(payloads{i},'Payload', 'DataSize'))/4,... 
      'int32','ieee-le');
end

%%% Read external signals
for i = 1:frame.nExternalSignals
   frame.ext(i).code = fread(fid,1,'int32','ieee-le');
   frame.ext(i).size = fread(fid,1,'int32','ieee-le');
   frame.ext(i).payload = fread(fid,frame.ext(i).size/4 - 2,'int32','ieee-le');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = readChars(fid,ln)
if nargin == 1, ln = 100; end;
out(1:ln) = ' ';
for i = 1:ln
   a = uint8(fread(fid,1,'uchar'));
   b = uint8(fread(fid,1,'uchar'));
   out(i) = char(bitor(bitshift(a,8),bitand(b,255)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = readUTF(fid)
ln = fread(fid,1,'int16');
out(1:ln) = ' ';
count = 1;
while count <= ln
    a  = uint8(fread(fid,1,'bit8'));
    if a < 128
        out(count) = char(a);
    elseif a < 224
        b = uint8(fread(fid,1,'bit8'));
        out(count) = char(bitor(shiftdim(bitand(a,31),6), bitand(b,63)));
    elseif a < 239
        b = fread(fid,1,'bit8');
        c = fread(fid,1,'bit8');
        out(count) = char(bitor(shiftdim(bitand(a,15),12), ...
                        bitor(shiftdim(bitand(b,63),6), ...
                           bitand(c,63))));
    end
    count = count + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = framevarsV3
% all the stuff in the frame
var = {
   'startIdent'               'int32' 'ieee-le'
   'endIdent'                 'int32'  'ieee-le'
   'type'                     'int32' 'ieee-le'
   'size'                     'int32' 'ieee-le'
   'timestamp'                'int32' 'ieee-le'
   'status'                   'int32'  'ieee-le'
   'errorCode'                'int32'  'ieee-le'
   'ecgDataSize'              'int32'  'ieee-le'
   'stretchDataSize'          'int32'  'ieee-le'
   'bendDataSize'             'int32'  'ieee-le'
   'positionDataSize'         'int32'  'ieee-le'
   'scanningPatternCodeDataSize' 'int32'  'ieee-le'
   'voltageInjectionDataSize' 'int32'  'ieee-le'
   'iqDataSize'               'int32'  'ieee-le'
   'nExternalSignals'         'int32'  'ieee-le'              
   };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = filevarsV3
% name and precision of the variables stored in the header
var = {
   'initialTimestamp'   'int64'
   'finalTimestamp'     'int64'
   'numberOfFrames'      'int32'
   'dataFileName'           'chars100'
   'dataFileCondition'      'chars300'
   'dataFileComments'       'chars600'
   'tic.imageRate'          'float32'
   'tic.injectionCurrent'   'float32'
   'tic.currentFreq'        'float32'
   'tic.switchSettlingTime' 'float32'
   'tic.injectionPattern'   'int32'
   'tic.num_elec'           'int32'
   'sbc.nco_freq'           'int32'
   'sbc.dac_gain'           'int32'
   'sbc.dac_fs'             'int32'
   'sbc.idta'               'charsNE'
   'sbc.mdta'               'charsNE'
   'sbc.tickl'              'int32'
   'sbc.tmckl'              'int32'
   'sbc.nppw'               'int32'
   'sbc.pga0_g'             'int32'
   'sbc.pga1_m'             'int32'
   'sbc.pga1_g'             'int32'
   'sbc.pga1_o'             'int32'
   'sbc.nign'               'int32'
   'sbc.tign'               'int32'
   'sbc.nsmp'               'int32'
   'sbc.tsmp'               'int32'
   'sbc.tr'                 'int32'
   'sbc.ni0'                'int32'
   'sbc.nq0'                'int32'
   'sbc.ni1'                'int32'
   'sbc.nq1'                'int32'
   'sbc.nraw'               'int32'
   'sbc.raw_m'              'int32'
   'sbc.endian'             'int32'
   };


