function [ TT  ] = ScouseTom_TrigProcess( Trigger,HDR )
%ScouseTom_TrigProcess Process trigger channels - Code to reject
%artefactual triggers and remove incomplete injections
%   Inputs:
%  Trigger  - struct from ScouseTom_TrigReadChn
% HDR - output from sopen


% This is separate from TrigReadChn as I will use the output of this but I
% can be bothered to argue with Kirill about how he does the EPs

% to do - clean fucked up injections switches


%% find only correct type of events in injection trigger channel

%sample rate
Fs=HDR.SampleRate;


switch HDR.TYPE
    case 'BDF' % biosemi file
        
        N_samples=HDR.NRec*Fs;
        
        
    case 'BrainVision'
        N_samples=HDR.SPR;
        
    otherwise
        error('Bad HDR');
end

%%
%define max period of INDENTIFICATION pulses at start of file, these are
%1000us apart
maxIDpulsemicros=2000; %max period of ID pulses to consider - this rejects all "real" pulses from start/stops/switching/stim
maxIDperiod = fix((maxIDpulsemicros*10^-6*Fs)); %rounded to nearest sample


%% Find the correct channels in data

Start_chn=find(strcmp(Trigger.Type, 'Start'));
Stop_chn=find(strcmp(Trigger.Type, 'Stop'));
Freq_chn=find(strcmp(Trigger.Type, 'Freq'));
Stim_chn=find(strcmp(Trigger.Type, 'Stim'));
Switch_chn=find(strcmp(Trigger.Type, 'Switch'));

if ~isempty(Start_chn)
    InjectionStarts=Trigger.RisingEdges{Start_chn};
else
    InjectionStarts=[];
end

if ~isempty(Stop_chn)
    InjectionStops=Trigger.RisingEdges{Stop_chn};
else
    InjectionStops=[];
end

if ~isempty(Switch_chn)
    Switches=Trigger.RisingEdges{Switch_chn};
else
    Switches=[];
end

if ~isempty(Stim_chn)
    Stims=Trigger.RisingEdges{Stim_chn};
else
    Stims=[];
end

if ~isempty(Freq_chn)
    Freqs=Trigger.RisingEdges{Freq_chn};
else
    Freqs=[];
end
%% check that all channels were read




%% Check is starts are missing

if isempty(InjectionStarts)
    %if there is no injection then add one at the very start of the file
    InjectionStarts(1)=1;
    disp('No Start Injection Found - Adding fake one at start of file');
end

%% Check Start codes

%Normal injections start with a single start pulse, contact checks start
%with two pulses

BelowThres = (diff(InjectionStarts) < maxIDperiod);
%use bwlabel to find connections in array
[S, NN]=bwlabel(BelowThres);

NumContact=NN; %number of contact starts in file
ContactStarts=InjectionStarts(S>0); %the contact start time is when the first pulse happens

%remove the start pulses which refer to the contact checks from the normal
%injection start vector

ContactStartsIdx=find(S > 0); % first pulse of contact ID
rem_idx=ContactStartsIdx+1; %its the SECOND pulse we want to remove
InjectionStarts(rem_idx)=[]; %get rid of them!

%these IDXs are then used after the injections have been segmented
ContactStartsIdx=find(ismember(InjectionStarts,ContactStarts));
InjectionStartsIdx=find(~ismember(InjectionStarts,ContactStarts));

NumInj=length(InjectionStartsIdx); %number of injections after removing contact starts

TotInj=length(InjectionStarts);

%% Output to user

fprintf('%d Injection starts and %d Contact starts found\n',NumInj,NumContact);

%% Process Each Injection

%for now treat each injection the same

InjectionSwitchesIN=cell(1,TotInj);
FreqChanges=InjectionSwitchesIN;
Stimulations=InjectionSwitchesIN;
ProtocolCompleteFlags=InjectionSwitchesIN;
FreqOrder=InjectionSwitchesIN;
FreqStarts=InjectionSwitchesIN;
FreqStops=InjectionSwitchesIN;


%% loop through each injection start and process the switches inside
for iInj=1:TotInj
    curStart=InjectionStarts(iInj); %current start
    curEnd= InjectionStops(find(InjectionStops > InjectionStarts(iInj),1,'First')); %find first stop after current start
    
    if isempty(curEnd) %if we didnt find any then fake one at end of file
        curEnd=N_samples -1;
        disp('No Stop Injection Found - Adding fake one at end of file');
        InjectionStops(iInj)=curEnd;
    end
    
    %find the indicators which belong to this injection
    InjectionSwitchesIN{iInj}= Switches (Switches >= curStart & Switches < curEnd);
    FreqChanges{iInj}= Freqs (Freqs >= curStart & Freqs < curEnd);
    Stimulations{iInj}= Stims (Stims >= curStart & Stims < curEnd);
    
    %% Process stuff
    %find the protocol complete flags
    [InjectionSwitchesIN{iInj},ProtocolCompleteFlags{iInj}]=findcompleteflags(InjectionSwitchesIN{iInj},maxIDperiod);
    
    %find the frequency order from freq pulses
    [FreqChanges{iInj},FreqOrder{iInj},FreqStarts{iInj}]=findfreqorder(FreqChanges{iInj},maxIDperiod);
    %arrange separate into freq starts and stops in matrix InjSwitches x
    %Freq - this makes processing easier
    [FreqStarts{iInj}, FreqStops{iInj},FreqOrder{iInj}]=reshapefreqtriggers(FreqChanges{iInj},FreqOrder{iInj},FreqStarts{iInj},InjectionSwitchesIN{iInj},InjectionStops(iInj),N_samples);
    
    % here is where you would do stim
    
    %% Put into more sensible form
    
    Nfreq=size(FreqStarts{1},2);
    
    
    InjectionSwitches=cell(TotInj,Nfreq);
    
    
    
    if isempty(FreqChanges{iInj})
        %if single frequency, then the windows come directly from
        %InjectionSwitches
        InjectionSwitches{iInj}=[InjectionSwitchesIN{iInj}, [InjectionSwitchesIN{iInj}(2:end); InjectionStops(iInj)]];
        
    else
        %if multifreq mode then we are interested in the Freq Starts and
        %Stops *not* the Injection Switches
        
        for iFreq=1:Nfreq
            
            %find the relevant index for the start and stops for this freq
            
            fqord_idx=FreqOrder{iInj} == iFreq;
            
            curSwitches=[sort(FreqStarts{iInj}(fqord_idx)), sort(FreqStops{iInj}(fqord_idx))];
            
            InjectionSwitches{iInj,iFreq}=curSwitches(all(~isnan(curSwitches),2),:);
        end
        
    end
    
end


%% Separate contact checks
% Make separate cell arrays for contact checks (so as not to confuse further processing)

Contact.InjectionSwitches=InjectionSwitches(ContactStartsIdx);
InjectionSwitchesIN(ContactStartsIdx)=[];
InjectionSwitches(ContactStartsIdx)=[];
Contact.FreqChanges=FreqChanges(ContactStartsIdx);
FreqChanges(ContactStartsIdx)=[];
Contact.Stimulations=Stimulations(ContactStartsIdx);
Stimulations(ContactStartsIdx)=[];
Contact.InjectionStarts=InjectionStarts(ContactStartsIdx);
InjectionStarts(ContactStartsIdx)=[];
Contact.InjectionStops=InjectionStops(ContactStartsIdx);




%% output dat shit

TT.InjectionSwitches =InjectionSwitches;
%TT.FreqChanges=FreqChanges;
TT.Stimulations=Stimulations;
TT.InjectionStops=InjectionStops;
TT.InjectionStarts=InjectionStarts;
TT.Trigger=Trigger; % store the trigger variable too
TT.Contact=Contact;
TT.ProtocolCompleteFlags=ProtocolCompleteFlags;
TT.FreqOrder=FreqOrder;
TT.FreqStarts=FreqStarts;
TT.FreqStops=FreqStops;

%we might want more detailed output here
disp('Triggers Processed OK');

end

function [InjSwitchOut,ProtCompFlag]=findcompleteflags(InjSwitchIn,maxIDperiod)

%finds the protocol complete flags - the double pulses on the switch
%channel

InjSwitchOut=InjSwitchIn;
BelowThres = (diff(InjSwitchIn) < maxIDperiod);
%use bwlabel to find connections in array
[S, NN]=bwlabel(BelowThres);

%%
ProtCompFlag=find(S > 0); % first pulse of contact ID
rem_idx=ProtCompFlag+1; %its the SECOND pulse we want to remove
InjSwitchOut(rem_idx)=[]; %get rid of them!

end


function [FreqChangesOut,FreqOrder,FreqStarts]=findfreqorder(FreqChangesIn,maxIDperiod)

%finds the freq order from the pulses on the frequency channel - these
%sequences of pulses only happen on the *start* of injection, so we can
%find the starts from these too



FreqChangesOut=FreqChangesIn;
BelowThres = (diff(FreqChangesIn) < maxIDperiod);
%use bwlabel to find connections in array
[S, NN]=bwlabel(BelowThres);

%% get freq order

%freq order is equal to the number of pulses in each set
FreqOrder=histc(S,1:NN);

FreqStarts=nan(NN,1);
FreqStartsIdx=FreqStarts;

for iFq=1:NN
    FreqStarts(iFq,1)=FreqChangesIn(find(S == iFq,1));
    FreqStartsIdx(iFq,1)=(find(S == iFq,1));
end

%% remove extra pulses so only those relating to start and stop in file
rem_idx=find(S > 0)+1; % any pulse found is bad,so shift up by one to get same idx as FreqCHangesIn
FreqChangesOut(rem_idx)=[]; %get rid of them!


end


function [FreqStarts, FreqStops,FreqOrderOut]=reshapefreqtriggers(FreqChanges,FreqOrderIn,FreqStartsIn,InjectionSwitches,InjectionStop,N_samples)

%find the corresponding freqstop pulse for each freq start, if one is
%missing due to the file being truncated then add a fake one at the end




%% find start and stops

FreqStopsIn=nan(size(FreqStartsIn));

for iFq=1:length(FreqStartsIn)
    
    curFstart=FreqStartsIn(iFq);
    curFstop=FreqChanges(find(FreqChanges > curFstart,1,'First'));
    
    %make fake stop if it is missing
    if isempty(curFstop)
        curFstop=InjectionStop-2;
    end
    
    FreqStopsIn(iFq,1)=curFstop;
    
end
%% arrange into each injection switch

N_freq=max(FreqOrderIn);
N_Sw=length(InjectionSwitches);

% nan pad as the last one might not be the correct size as it might have
% been truncated
FreqStarts=nan(N_Sw,N_freq);
FreqStops=FreqStarts;
FreqOrderOut=FreqStarts;


for iSw=1:N_Sw
    
    %find the index of frequency starts within the current injection switch
    
    if iSw < N_Sw
        cur_idx=find(((FreqStartsIn > InjectionSwitches(iSw) & FreqStartsIn <= InjectionSwitches(iSw+1))));
    else % there is no injection switch for the final stop - so use injection stop
        cur_idx=find(((FreqStartsIn > InjectionSwitches(iSw) & FreqStartsIn <= InjectionStop)));
    end
    %I found that due to the speed of things on the Due, the Switch
    %Indication could actually occur *before* the frequency stop indicator (though the ISR),
    %despite coming ~50 lines of code afterwards - which means just finding
    %the stop flags within the injection would sometimes not work. this
    %took *way* too long to figure out.
       
    FreqStarts(iSw,1:length(cur_idx))=FreqStartsIn(cur_idx);
    FreqStops(iSw,1:length(cur_idx))=FreqStopsIn(cur_idx);
    FreqOrderOut(iSw,1:length(cur_idx))=FreqOrderIn(cur_idx);
    
    
    
end




end



















































