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




%% Separate into injections

if isempty(InjectionStarts)
    %if there is no injection then add one at the very start of the file
    InjectionStarts(1)=1;
    disp('No Start Injection Found - Adding fake one at start of file');
end

NumInj=length(InjectionStarts);

InjectionSwitches=cell(1,NumInj);
FreqChanges=InjectionSwitches;
Stimulations=InjectionSwitches;

for iInj=1:NumInj
    curStart=InjectionStarts(iInj); %current start
    curEnd= InjectionStops(find(InjectionStops > InjectionStarts,1,'First')); %find first stop after current start
    
    if isempty(curEnd) %if we didnt find any then fake one at end of file
        curEnd=N_samples -1;
        disp('No Stop Injection Found - Adding fake one at end of file');
    end
    
    %find the indicators which belong to this injection
    InjectionSwitches{iInj}= Switches (Switches >= curStart & Switches < curEnd);
    FreqChanges{iInj}= Freqs (Freqs >= curStart & Freqs < curEnd);
    Stimulations{iInj}= Stims (Stims >= curStart & Stims < curEnd);
    
    %Clean the injections somehow....
    
    
end

%% output dat shit
disp('Injection channels read just fine');

TT.InjectionSwitches =InjectionSwitches;
TT.FreqChanges=FreqChanges;
TT.Stimulations=Stimulations;
TT.InjectionStops=InjectionStops;
TT.InjectionStarts=InjectionStarts;
TT.Trigger=Trigger; % store the trigger variable too


end
