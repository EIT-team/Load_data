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
%number of samples in file
N_samples=HDR.NRec*Fs;

%look for events only on the current channel - chn 0

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
rem_idx=ContactStartsIdx+1; %its the SECOND pulse we want to remove pulse
InjectionStarts(rem_idx)=[]; %get rid of them!

%these IDXs are then used after the injections have been segmented
ContactStartsIdx=find(ismember(InjectionStarts,ContactStarts));
InjectionStartsIdx=find(~ismember(InjectionStarts,ContactStarts));

NumInj=length(InjectionStartsIdx); %number of injections after removing contact starts

TotInj=length(InjectionStarts);

%% Output to user

fprintf('%d Injection starts and %d Contact starts found\n',NumInj,NumContact);


%% Process Each Injections

%for now treat each injection the same

InjectionSwitches=cell(1,TotInj);
FreqChanges=InjectionSwitches;
Stimulations=InjectionSwitches;

for iInj=1:TotInj
    curStart=InjectionStarts(iInj); %current start
    curEnd= InjectionStops(find(InjectionStops > InjectionStarts,1,'First')); %find first stop after current start
    
    if isempty(curEnd) %if we didnt find any then fake one at end of file
        curEnd=N_samples -1;
        disp('No Stop Injection Found - Adding fake one at end of file');
        InjectionStops(end)=curEnd;
    end
    
    %find the indicators which belong to this injection
    InjectionSwitches{iInj}= Switches (Switches >= curStart & Switches < curEnd);
    FreqChanges{iInj}= Freqs (Freqs >= curStart & Freqs < curEnd);
    Stimulations{iInj}= Stims (Stims >= curStart & Stims < curEnd);
    
    %Clean the injections somehow....
    
end


%% Separate contact checks
% Make separate cell arrays for contact checks (so as not to confuse further processing)

Contact.InjectionSwitches=InjectionSwitches(ContactStartsIdx);
InjectionSwitches(ContactStartsIdx)=[];
Contact.FreqChanges=FreqChanges(ContactStartsIdx);
FreqChanges(ContactStartsIdx)=[];
Contact.Stimulations=Stimulations(ContactStartsIdx);
Stimulations(ContactStartsIdx)=[];
Contact.InjectionStarts=InjectionStarts(ContactStartsIdx);
InjectionStarts(ContactStartsIdx)=[];
Contact.InjectionStops=InjectionStops(ContactStartsIdx);




%% output dat shit
disp('Triggers Processed OK');

TT.InjectionSwitches =InjectionSwitches;
TT.FreqChanges=FreqChanges;
TT.Stimulations=Stimulations;
TT.InjectionStops=InjectionStops;
TT.InjectionStarts=InjectionStarts;
TT.Trigger=Trigger; % store the trigger variable too
TT.Contact=Contact;


end
