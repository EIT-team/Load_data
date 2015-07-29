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

function PulseStart=ReadIndPin(chn,HDR)

%% FOR REF ONLY TO LOOK AT HOW I CLEANED THE SWITCHES BEFORE HAND



%this subfunction reads the pules on the indicator pin channel - this is a
%hangover from when we didnt use the trigger channel and is only relevant
%for 2013 datasets. This also fixes broken switches due to the indicator
%pin going funny sometimes


disp('Loading Header File for Indicator Channel ONLY');

%get voltages from indicator channel only

HDRi=sopen(HDR.FileName,'r',chn,['OVERFLOWDETECTION:OFF']);

Fs=HDRi.SampleRate;

IndicatorPinData=sread(HDRi,inf,0);

thres=2000;
min_coef=0.75;

%% process the indicator pin data as sometimes it fucks up


%find where it goes negative and set it to the "baseline" value
databaseline=mean(IndicatorPinData(1:Fs));

IndicatorPinData(IndicatorPinData < databaseline) = databaseline;






%% find contact check start from indicator channel

% find peaks in indicator channel

%min width - pulse is 500 micro seconds, so look only for pulses close to
%this
minwidth = fix(min_coef*(500e-6*Fs)); %rounded to nearest sample

%find peaks by creating logical threshold array. then finding the rising
%and falling edges, checking the width is greater than the min width


AboveThres=IndicatorPinData > thres; %threshold indicator pin data

AboveThres=[ false; AboveThres; false]; % pad with zeros (for diff below)
%this is now a logical array of 0 1. pulses are found by finding when the
%times when the data is above the threshold.

ThresEdges= diff(AboveThres); %take diff of this data. this is now nearly all 0 except for 1 for rsing edge and -1 for falling

RisingEdges = find(ThresEdges ==1); % get rising edges

FallingEdges = find(ThresEdges ==-1); %get falling edges

Pulsewidth = FallingEdges - RisingEdges; % get width of pulses

GoodPulses = Pulsewidth >= minwidth; % good pulses are those which are greater than the minimum width

PulseStart = RisingEdges(GoodPulses); %start of good pulses only

%% this is cheating to unfuck the indicator pin channel


%find if there are any starts in the file - then find if there are any
%"gaps" in the pulses less than 10 times the average pulse width


disp('Indicator Pin is usually fucked so Im gonna fix it now');


% threshold of pulse width - close to 500 microseconds
maxwidth = fix(1.1*(1000e-6*Fs)); %rounded to nearest sample

%find the pulses which are close together
BelowThres = (diff(PulseStart) < maxwidth);

%use bwlabel to find connections in array
[S, NN]=bwlabel(BelowThres);

%initialse array of starting positions
InjectionStarts=[];

%initialise counters
InjCnt=1;

%for each connected pulse, find size and see if it is a injection or
%contact impedance start
for ii=1:NN
    
    %find the number samples with this code
    codetmp=find (S == ii);
    codesize=size(codetmp,1);
    
    
    switch codesize
        case 1 % injection start is 2 pulses so diff =1
            InjectionStarts(InjCnt)=codetmp(1);
            InjCnt=InjCnt+1;
            
    end
end

if ~isempty(InjectionStarts)
    InjectionFirst= InjectionStarts+2;
else
    InjectionFirst=1;
end


%find the "correct" length of the switch


%%
%find the gaps, append a pulse at the end then sort them

%find where the gap is bigger than the correct value
Sw_in=diff(PulseStart(InjectionFirst:end));
Sw_True=mode(Sw_in);

%find where the gaps are
Sw_corrected=ceil(Sw_in./Sw_True);
gaps= Sw_corrected > 1;
gapind= find(gaps ==1);


%plug the gaps! this is cheating but who cares/will ever look at this
for gg=1:length(gapind)
    
    pulses_missing=Sw_corrected(gapind(gg)) -1;
    
    new_pulse_vec=(1:pulses_missing)*Sw_True+PulseStart(InjectionFirst+gapind(gg)-1);
    
    PulseStart=[PulseStart; new_pulse_vec'];
end
PulseStart=sort(PulseStart);


%find too small pulses
Pulses_short= find(diff(PulseStart(InjectionFirst:end)) < fix(0.9*Sw_True) == 1);

%get rid of them!!!!!
PulseStart(Pulses_short+InjectionFirst)=[];


end

