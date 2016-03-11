function [ Trigger ] = ScouseTom_TrigReadChn( HDR,SkipIDCodes )
%SCOUSETOM_READTRIGCHN Identifies events on trigger channels, and identify
%which is which according to the ID codes at the start of the file
%   Detailed explanation goes here

%% Get trigger channel input according to which file type it is
switch HDR.TYPE
    case 'BDF' % biosemi file
        trignum=8;
        [ StatusChns,TrigPos ] = ScouseTom_getbdftrig( HDR,trignum );
        % Other useful info from HDR
        
    case 'BrainVision'
        [ StatusChns,TrigPos ] = ScouseTom_geteegtrig( HDR );
        fname=HDR.FILE.Name;
        trignum=size(StatusChns,2);
    otherwise
        error('Bad HDR');
end

Fs=HDR.SampleRate;

%% Define Variables

if exist('SkipIDCodes','var') ==0
    SkipIDCodes=0;
end

%number of trigger channels
trignum=8; % 8 for BioSemi and ActiChamp

%threshold value for finding edges
thres=0.5; % can be so high as data *MUST* be logical 1 and 0

%define min width of pulses - to reject spurious noisey ones
minpulsemicros=150; % anything less than 150 does not count
minwidth = floor(((minpulsemicros*10^-6)*Fs)); %rounded to nearest sample

%define max period of INDENTIFICATION pulses at start of file, these are
%1000us apart
maxIDpulsemicros=2000; %max period of ID pulses to consider - this rejects all "real" pulses from start/stops/switching/stim
maxIDperiod = fix((maxIDpulsemicros*10^-6*Fs)); %rounded to nearest sample


%% Channel ID Codes

%these are the number of GAPS BETWEEN PULSES *NOT* Pulses
ID_Codes.Name(1)={'Stim'};
ID_Codes.Num(1)=1; % Stim is 2 pulses so diff =1
ID_Codes.Name(2)={'Start'};
ID_Codes.Num(2)=2; % injection start is 3 pulses so diff =2
ID_Codes.Name(3)={'Stop'};
ID_Codes.Num(3)=4;
ID_Codes.Name(4)={'Switch'};
ID_Codes.Num(4)=3;
ID_Codes.Name(5)={'Freq'};
ID_Codes.Num(5)=5;



ID_Codes.DefaultOrder=[3,1,2,


%there may be others here - system has 3 spare channles EX_1 2 and 3 on
%arduino. and Kirills physchotool box stuff will also go here

%% CHECK HDR IS OK here



%% FIND EDGES IN EACH CHANNEL

%find peaks by creating logical threshold array. then finding the rising
%and falling edges, checking the width is greater than the min width

AboveThres=StatusChns > thres; %threshold indicator pin data

AboveThres=[ AboveThres(1,:); AboveThres; AboveThres(end,:);]; % pad array (for diff below)
%this is now a logical array of 0 1. pulses are found by finding when the
%times when the data is above the threshold.

ThresEdges= diff(AboveThres); %take diff of this data. this is now nearly all 0 except for 1 for rsing edge and -1 for falling

[RisingEdgesPosIdx,RisingEdgesChn] = find(ThresEdges ==1); % get rising edges

[FallingEdgesPosIdx,FallingEdgesChn] = find(ThresEdges ==-1); %get falling edges

RisingEdges=TrigPos(RisingEdgesPosIdx);
FallingEdges=TrigPos(FallingEdgesPosIdx);


%% read triggers in each channel and reject orphaned ones or too short ones

%find any orphaned falling edges at start of file which can happen if the
%biosemi cries a little bit or if the recording was started *before*
%arduino turns on. As pins are held HIGH by biosemi. And find orphaned
%rising too, which can occur if the file stops early

for iChn=1:trignum
    %take only rising and falling belonging to this channel
    curRising=RisingEdges(RisingEdgesChn == iChn);
    curFalling=FallingEdges(FallingEdgesChn == iChn);
    
    if (~isempty(curRising) && ~isempty(curFalling)) % if R and F edges found
        
        if curFalling(1) < curRising(1) %if first falling happening before first rising
            curFalling(1)=[];
        end
        
        if curRising(end) > curFalling(end)
            curRising(end) =[];
        end
        
    else %if both R and F dont exist then ignore all edges in this channel
        
        curRising=[];
        curFalling=[];
        
    end
    
    %now we have only actual pulses in the channel, remove short/weird
    %pulses
    
    Pulsewidth = curFalling - curRising; % get width of pulses
    GoodPulses = Pulsewidth >= minwidth; % good pulses are those which are greater than the minimum width
    Trigger.RisingEdges(iChn)={curRising(GoodPulses)};
    Trigger.FallingEdges(iChn)={curFalling(GoodPulses)};
    
end


%% NEXT IDENTIFY CHANNELS BY READING THE LITTLE COMMAND ONES TO START WITH


if ~SkipIDCodes

%counter for unknown trigger channels
ChnUnknown=0;

for iChn=1:trignum
    
    PulseStart=Trigger.RisingEdges{iChn};
    %find the pulses which are close together
    BelowThres = (diff(PulseStart) < maxIDperiod);
    %use bwlabel to find connections in array
    [S, NN]=bwlabel(BelowThres);
    
    %assume only 1 ID and this happens first
    
    codetmp=find (S == 1);
    codesize=size(codetmp,1);
    
    %find the ID code this belong to
    chnID=find(ID_Codes.Num == codesize);
    
    if ~isempty(chnID)
        
        %assign this trigger a name
        Trigger.Type(iChn)=ID_Codes.Name(chnID);
        Trigger.ID_Code(iChn)=codesize;
        
        %delete the references to the rising and falling edges for the ID codes
        rem_idx=[true; S];
        rem_idx=find(rem_idx ==1);
        
        Trigger.ID_Rising(iChn)={Trigger.RisingEdges{iChn}(rem_idx)};
        Trigger.ID_Falling(iChn)={Trigger.FallingEdges{iChn}(rem_idx)};
        
        Trigger.FallingEdges{iChn}(rem_idx)=[];
        Trigger.RisingEdges{iChn}(rem_idx)=[];
        
    else
        
        if ~isempty(PulseStart)
            Trigger.Type(iChn)={['Unknown_' num2str(ChnUnknown)]}; %create name
            ChnUnknown=ChnUnknown+1;
            Trigger.ID_Code(iChn)=0;
        else
            Trigger.ID_Code(iChn)=nan;
            Trigger.Type(iChn)={''};
        end
    end
    
end

else
    fprintf(2,'SKIPPING ID CODE CHECK - ASSUMING EVERYTHING WIRED CORRECTLY!');
    Trigger.ID_Code=[3,1,2,4,5,nan,nan,nan];
    
end


%% output little bit of info about what was found

%find channels that actually had pulses in them
GoodChn=cellfun(@(x) ~isempty(x),Trigger.RisingEdges);

UnknownChn=cellfun(@(x) ~isempty(x),strfind(Trigger.Type,'Unknown'));

%get rid of Unkown channels in "good channels"
GoodChn=GoodChn ~= UnknownChn;

NumUnknownChn=sum(UnknownChn);
NumGoodChn=sum(GoodChn);

fprintf('Found %d trig chn: ',NumGoodChn);
fprintf('%s, ',Trigger.Type{GoodChn});
fprintf('and %d unknown channel(s),',NumUnknownChn);
fprintf(' so thats good.\n');


end

