function [ Trigger ] = ScouseTom_TrigReadChn( HDR,IdentifyChannels,SkipIDCodes,TimeToIgnore)
% [ Trigger ] = ScouseTom_TrigReadChn( HDR,IdentifyChannels,SkipIDCodes,TimeToIgnore)
% SCOUSETOM_READTRIGCHN Identifies events on trigger channels, and identify
% which is which according to the ID codes at the start of the file
%
% This code has 3 jobs: ID trigger channels, remove too short pulses, sort
% the pulses into the different channels for later processing. This is
% necessary as the BioSemi and ActiChamp give different types of infomation
% about the digital channel.
%
%
% Inputs [default]:
% HDR - HDR from ScouseTom_getHDR

% Identify channels [def 0] -  The ScouseTom system writes ID codes at the
% start of an injection to ID which channels are which.
% This is because at the time they were still connected manually,
% and the pins would get swapped. This is not needed if you are using the
% proper shield, so this is just skipped by default.

% SkipIDCodes [def 0] - Ignore checking for ID codes at all, and leave all
% channels unlabelled. This is useful if you are not using the normal
% triggers, like with the robot arm in the tank. 

% TimeToIgnore[def 0] - Time in seconds which will be removed from
% processing. Useful if you have some unusual triggers, or the arduino
% reset during recording. 

% Output:
% Trigger - Structure used in subsequent processing steps
% ScouseTom_TrigProcess. Containing Cells of each trigger type - Start Stop
% Switch Stim etc.

%% Get trigger channel input according to which file type it is
switch HDR.TYPE
    case 'BDF' % biosemi file
        trignum=8;
        [ StatusChns,TrigPos ] = ScouseTom_getbdftrig( HDR,trignum );
        % Other useful info from HDR
        
    case 'BrainVision'
        [ StatusChns,TrigPos ] = ScouseTom_geteegtrig( HDR );
        
    otherwise
        error('Bad HDR');
end

Fs=HDR.SampleRate;

%% Define Variables
% by default we do not identify the channels, and assume the numbering was
% correct
if exist('IdentifyChannels','var') ==0  || isempty(IdentifyChannels)
    IdentifyChannels=0;
end
% by deafult we do not skip finding the ID code blocks, as they are
% normally there in the recording. But we do want to remove them for later
% stages, as they do not represent "useful" infomation
if exist('SkipIDCodes','var') ==0
    SkipIDCodes=0;
end

%number of trigger channels
trignum=8; % 8 for BioSemi and ActiChamp

%threshold value for finding edges
thres=0.5; % can be so high as data *MUST* be logical 1 and 0

%define min width of pulses - to reject spurious noisey ones
% minpulsemicros=150; % anything less than 150 does not count
minpulsemicros = 1; % Set this to one as for some reason some actichamp ID codes are only 1 sometimes!
minwidth = max([floor(((minpulsemicros*10^-6)*Fs)) 1]); %rounded to nearest sample at least 1

%define max period of INDENTIFICATION pulses at start of file, these are
%1000us apart
maxIDpulsemicros=2000; %max period of ID pulses to consider - this rejects all "real" pulses from start/stops/switching/stim
maxIDperiod = ceil(((maxIDpulsemicros*10^-6)*Fs)); %rounded to nearest sample


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
ID_Codes.Name(6)={'Compliance'};
ID_Codes.Num(6)=6;
ID_Codes.Name(7)={'DummyIgnoreThisOne'};
ID_Codes.Num(7)=8;
ID_Codes.Name(8)={'Unknown'};
ID_Codes.Num(8)=-1;

ID_Codes.DefaultOrder=[4,1,2,3,5,6,8,7];
ID_Codes.DefaultID=ID_Codes.Num([ID_Codes.DefaultOrder]);
ID_Codes.DefaultName=ID_Codes.Name([ID_Codes.DefaultOrder]);

ID_Codes.DefaultOrder(end+1:trignum)=nan;
ID_Codes.DefaultID(end+1:trignum)=nan;
ID_Codes.DefaultName(end+1:trignum)={''};

%% CHECK HDR IS OK here



%% DELETE ONES WE DONT WANT

if exist('TimeToIgnore','var')
    
    rem_idx = (TrigPos/Fs) < TimeToIgnore;
    
    TrigPos(rem_idx) =[];
    StatusChns(rem_idx,:) =[];
    
end


%% FIND EDGES IN EACH CHANNEL

%only biosemi has proper record of rising AND falling edges, so we need to
%process them slightly differently. For the diff to find Thresedges, we can
%pad with the edge values for biosemi. But for the ActiChamp we need to pad
%with zeros to ensure we have a rising edge at start


%find peaks by creating logical threshold array. then finding the rising
%and falling edges, checking the width is greater than the min width

AboveThres=StatusChns > thres; %threshold indicator pin data

% AboveThres=[ AboveThres(1,:); AboveThres; AboveThres(end,:);]; % pad array (for diff below)
switch HDR.TYPE
    case 'BDF' % biosemi file
        AboveThres=[ AboveThres(1,:); AboveThres; AboveThres(end,:);]; % pad array (for diff below)
        
    case 'BrainVision'
        AboveThres=[ zeros(size(AboveThres(1,:))); AboveThres; AboveThres(end,:);]; % pad array (for diff below)
        
    otherwise
        error('Bad HDR');
end


%this is now a logical array of 0 1. pulses are found by finding when the
%times when the data is above the threshold.

ThresEdges= diff(AboveThres); %take diff of this data. this is now nearly all 0 except for 1 for rsing edge and -1 for falling

[RisingEdgesPosIdx,RisingEdgesChn] = find(ThresEdges ==1); % get rising edges

[FallingEdgesPosIdx,FallingEdgesChn] = find(ThresEdges ==-1); %get falling edges

RisingEdges=TrigPos(RisingEdgesPosIdx);
FallingEdges=TrigPos(FallingEdgesPosIdx);


%% read triggers in each channel and reject orphaned ones or too short ones

% only BioSemi has rising AND falling edges, so we cannot do this with
% actichamp


%find any orphaned falling edges at start of file which can happen if the
%biosemi cries a little bit or if the recording was started *before*
%arduino turns on. As pins are held HIGH by biosemi. And find orphaned
%rising too, which can occur if the file stops early

for iChn=1:trignum
    %take only rising and falling belonging to this channel
    curRising=RisingEdges(RisingEdgesChn == iChn);
    curFalling=FallingEdges(FallingEdgesChn == iChn);
    
    %     figure;hold on;stairs(TrigPos,StatusChns(:,iChn));plot(curRising,[1],'o');plot(curFalling,[1],'x');hold off;title(['chn : ' num2str(iChn)]);
    
    
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

%counter for unknown trigger channels
ChnUnknown=0;

% check the ID codes at start of file to find which channels are which -
% but we can assume its wired in standard way in nearly all cases
if  IdentifyChannels
    
    
    for iChn=1:trignum
        
        %assume all ID codes must happen in the first second after first
        %trigger. this is to limit multiple pulses like freq changes or
        %protocol complete being read as id codes
        
        PulseStart=Trigger.RisingEdges{iChn} (Trigger.RisingEdges{iChn}< RisingEdges(1)+Fs);
        %         PulseStart
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
    fprintf(2,'SKIPPING ID CODE CHECK - ASSUMING EVERYTHING WIRED CORRECTLY!\n');
    Trigger.ID_Code=ID_Codes.DefaultID;
    Trigger.Type=ID_Codes.DefaultName;
    
end


%% Identify blocks of ID Codes
% The start of each injection creates an ID Code block. Search each of them
% in turn in case there is multiple EIT injections in a single recording.

if ~SkipIDCodes
    
    
    %% Using the Start channel find groups
    StartChn=find((cellfun(@(x) ~isempty(x),strfind(Trigger.Type,ID_Codes.Name{2}))));
    
    % This time take *all* pulses
    PulseStart=Trigger.RisingEdges{StartChn};
    
    %find the pulses which are close together
    BelowThres = (diff(PulseStart) < maxIDperiod);
    %use bwlabel to find connections in array
    [S, NN]=bwlabel(BelowThres);
    
    %check groups are the correct size
    for iGroup = 1:NN
        if size(find (S == iGroup),1) ~= Trigger.ID_Code(StartChn)
%             error('Start ID not correct on Start channel!');
            S(S==iGroup) =nan;
            NN=NN-1;
            S(S>iGroup) = S(S>iGroup)-1;
            S(isnan(S))=[];
            
        end
    end
    
    fprintf('Found %d Injection starts\n',NN);
    
    
    % find all ID codes around this window
    
    for iGroup = 1:NN
        % find where the ID codes start
        IDCodeStart = PulseStart(find(S ==iGroup,1));
        
        for iChn=1:trignum
            % find all rising edges in a window around this starting pulse
            rem_idx = find( Trigger.RisingEdges{iChn} > IDCodeStart - maxIDperiod*10 & Trigger.RisingEdges{iChn} < IDCodeStart + maxIDperiod*10);
            % store them as reference
            Trigger.ID_Rising(iChn)={Trigger.RisingEdges{iChn}(rem_idx)};
            Trigger.ID_Falling(iChn)={Trigger.FallingEdges{iChn}(rem_idx)};
            
            %remove them from the main array
            Trigger.FallingEdges{iChn}(rem_idx)=[];
            Trigger.RisingEdges{iChn}(rem_idx)=[];
            
            
        end
        
    end
    
else
    % the Trigger structure is not altered
end




%% Clear the dummy channel

%we have extra channel to force rising edges so the actichamp triggers make
%sense. But we dont want this data so clear the channel


%find the one with the correct ID code - This works for BDF files or if we
%have forced defaults above
if any(cellfun(@(x) ~isempty(x),strfind(Trigger.Type,ID_Codes.Name{7})))
    
    DummyChn=find((cellfun(@(x) ~isempty(x),strfind(Trigger.Type,ID_Codes.Name{7}))));
    
else
    %find it from the first 8 rising edges - This works for ActiChamp
    
    %find first event on chn other than the dummy one
    MainIDCodeStart=find(sum(StatusChns,2) >1,1);
    
    %find what channels have rising edges in them - this should *all* be the
    %dummy channel (usually 8)
    [~, startidchn]=find(StatusChns(1:MainIDCodeStart-1,:));
    
    %if any of the first pulses happen *at least* 75% on the same channel, then
    %that is the dummy channel
    DummyChn=find(histc(startidchn,1:trignum) > length(startidchn)*0.75);
end

if ~isempty(DummyChn)
    
    
    Trigger.Type(DummyChn)={''};
    Trigger.ID_Code(DummyChn)=nan;
    Trigger.RisingEdges(DummyChn)={[]};
    Trigger.FallingEdges(DummyChn)={[]};
    Trigger.ID_Rising(DummyChn)={[]};
    Trigger.ID_Falling(DummyChn)={[]};
    
    
end

%% Check if ok and output

%at the moment, if we dont find switch, start and stop *at least* then
if (any(Trigger.ID_Code == ID_Codes.Num(2)) ...
        && any(Trigger.ID_Code == ID_Codes.Num(3)) ...
        && any(Trigger.ID_Code == ID_Codes.Num(4)))
    
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
    
else
    
    if ~SkipIDCodes
        fprintf(2,'STARTING CODES WERE BROKEN! Trying again but assuming start of file missing, and forcing default channels\n');
        
        Trigger=ScouseTom_TrigReadChn(HDR,1);
    else
        fprintf(2,'TRIGGERS ARE MESSED UP! Giving up \n');
    end
    
    
    
end

