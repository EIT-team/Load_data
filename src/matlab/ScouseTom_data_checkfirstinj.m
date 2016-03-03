function [ StartInj ] = ScouseTom_data_checkfirstinj(HDR,InjectionSwitchesCell,Protocol )
%SCOSUETOM_DATA_CHECKFIRSTINJ Checks the first injection in the dataset to
%see if it is the correct protocol line, adjusts the processing if not.
%This is to account for
%   Detailed explanation goes here

%% check inputs are ok

%if the data is fucked or has big artefact then this can mess up as the RMS
%is borked



%% checking for the double prt pulse goes here



%% Load first second of data

Nfreq=size(InjectionSwitchesCell,2);
Fs=HDR.SampleRate;

StarInj=nan(Nfreq,1);

%find earliest switch
for iFreq=1:Nfreq
   
    FirstSample=min(InjectionSwitchesCell{iFreq}(1,:));
    StartSec=floor(FirstSample/Fs);
    StartSample=StartSec*Fs;
    
    %% Load a chunk of data
    
    %only load 2 seconds as we want *at most* 1 second for the estimation
    V=sread(HDR,2,StartSec);
    
    %% find which bit of the dataset to use
    
    %take either the first injection or the first second
    
    tmp=InjectionSwitchesCell{iFreq}(1,2)-InjectionSwitchesCell{iFreq}(1,1);
    if tmp > Fs
        tmp=Fs; %if the first switch is longer than a second, only take a second
    end
    
    
    tmpstart=InjectionSwitchesCell{iFreq}(1,1)-StartSample;
    tmpidx=tmpstart:tmpstart+tmp;
    
    
    
    %% estimate the injection pair
    
    Threshold=0.6; %coefficient for deciding good injection pair estimate
    
    %estimate the injection pairs from the two largest RMS values
    [StartInjEst, estimatebadness,RmsEst]=ScouseTom_data_EstInjPair(V(tmpidx,:),Threshold);
    
    %find desired first startinj - sort as we cant tell sources from sinks
    Protocol=sort(Protocol,2);
    
    %get the injection pair from the protocol
    StartInjProt=Protocol(1,:)';
    
    %find matching line in protocol
    start_poss=find(all([StartInjEst(1)==Protocol(:,1) StartInjEst(2)==Protocol(:,2)],2));
    
    if estimatebadness ==1
        fprintf(2,'WARNING! The Injection start estimate did not meet threshold\n');
    end
    
    
    if ~isempty(start_poss)
        fprintf('Starting injection pair for freq %d was found to be : %d\n',iFreq, start_poss);
        StartInj(iFreq)=start_poss;
    else
        fprintf(2,'NO MATCHING INJECTION PAIR FOUND! DID YOU LOAD CORRECT EXPSETUP?\nCheck bar plot and consider changing threshold');
        
        figure;
        hold on
        bar(RmsEst)
        plot([0 length(RmsEst)],repmat([max(RmsEst)*Threshold],1,2),'k-','Linewidth',2)
        xlabel('Channel');
        ylabel('V RMS in first injection');
        
        StartInj(iFreq)=-1;
        estimatebadness=1;
        
    end
    
end



end

