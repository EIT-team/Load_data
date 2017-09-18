function [ StartInj,EstimateWarning,RMSEst ] = ScouseTom_data_checkfirstinj(HDR,InjectionSwitchesCell,Protocol )
%[StartInj,EstimateWarning,RMSEst] = ScouseTom_data_checkfirstinj(HDR,InjectionSwitchesCell,Protocol )
%SCOUSETOM_DATA_CHECKFIRSTINJ Checks the first injection in the dataset and
%matches it to a line in the current injection protocol.
%
%
% Takes the RMS of the voltages given, finds the two largest channels, and
% matches these to an injection pair in the protocol given. If one is not
% found then shows the RMS and prompts user. If this is the case then
% normally something is broken in the data
%
%   Inputs:
%   HDR - from ScouseTom_getHDR
%   InjectionSwitchesCell - From the TT strucutre from
%       ScouseTom_TrigProcess. Normally given by TT.InjectionSwitches(TT.InjectionSwitches(1,:)
%   Protocol - Array of injection pairs usually ExpSetup.Protocol
%
%   Outputs:
%   StartInj(nFreq) - Line in protocol which each frequency begins at
%   EstimateWarning - Whether there was a clear injection pair. Threshold
%       is 3rd largest channel is 0.6 of 2nd.
%   RMSEst - the rms on all channels used in estimating inj pairs
%% check inputs are ok





%% checking for the double prt pulse goes here



%% Info from inputs

Nfreq=size(InjectionSwitchesCell,2);
Fs=HDR.SampleRate;

StartInj=nan(Nfreq,1);

%% Find start injection for each frequency

for iFreq=1:Nfreq
    %find earliest switch
    FirstSample=min(InjectionSwitchesCell{iFreq}(1,:));
    StartSec=floor(FirstSample/Fs); %findnearest second
    StartSample=StartSec*Fs; %corresponding sample
    
    %% Load a chunk of data
    
    %only load 2 seconds as we want *at most* 1 second for the estimation
    V=sread(HDR,2,StartSec);
    
    %% find which bit of the dataset to use
    
    %take either the first injection or the first second
    
    tmp=InjectionSwitchesCell{iFreq}(1,2)-InjectionSwitchesCell{iFreq}(1,1);
    if tmp > Fs
        tmp=Fs; %if the first switch is longer than a second, only take a second
    end
    
    %index of V we actually want
    tmpstart=InjectionSwitchesCell{iFreq}(1,1)-StartSample;
    tmpidx=tmpstart:tmpstart+tmp;
    
    
    %% estimate the injection pair
    
    Threshold=0.6; %coefficient for deciding good injection pair estimate
    
    %estimate the injection pairs from the two largest RMS values
    [StartInjEst, EstimateWarning,RMSEst]=ScouseTom_data_EstInjPair(V(tmpidx,:),Threshold);
    
    %find desired first startinj - sort as we cant tell sources from sinks
    Protocol=sort(Protocol,2);
    
    %get the injection pair from the protocol
    StartInjProt=Protocol(1,:)';
    
    %find matching line in protocol
    start_poss=find(all([StartInjEst(1)==Protocol(:,1) StartInjEst(2)==Protocol(:,2)],2));
    
    %warn if estimate was not matching threshold
    if EstimateWarning ==1
        fprintf(2,'WARNING! The Injection start estimate did not meet threshold\n');
    end
    
    %if one was found then hooray
    if ~isempty(start_poss)
        fprintf('Found start injection pair for freq %d: %d\n',iFreq, start_poss);
        StartInj(iFreq)=start_poss;
    else
        %if one was not found, then plot graph
        fprintf(2,'NO MATCHING INJECTION PAIR FOUND! DID YOU LOAD CORRECT EXPSETUP?\nCheck bar plot and consider changing threshold\n');
        
        figure;
        hold on
        bar(RMSEst)
        plot([0 length(RMSEst)],repmat([max(RMSEst)*Threshold],1,2),'k-','Linewidth',2)
        xlabel('Channel');
        ylabel('V RMS in first injection');
        title(['RMS in first inj. Freq ' num2str(iFreq) ', expected ' num2str(StartInjProt(1)) ' & ' num2str(StartInjProt(2)) ])
        
        StartInj(iFreq)=-1;
        EstimateWarning=1;
        
    end
    
end



end

