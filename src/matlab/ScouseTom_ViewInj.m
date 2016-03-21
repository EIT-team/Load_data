function [  ] = ScouseTom_ViewInj( HDR,ChntoView,Repeat,FreqNum,StartInj,TT,ExpSetup )
%SCOUSETOM_VIEWINJ Summary of this function goes here
%   Detailed explanation goes here

switch HDR.TYPE
    case 'BDF' % biosemi file
        %biosemi is stored as 1sec long records. so number of seconds is
        %Nrec
        Nsec=HDR.NRec;
        %biosemi file size is stored in bytes
        Fsize=HDR.FILE.size;
    case 'BrainVision'
        Nsec=ceil(HDRac.SPR/HDRac.SampleRate);
        %actichamp HDR points to header .vhdr file, so filesize is wrong
        eegfile=dir(fullfile(HDR.FILE.Path,[HDR.FILE.Name '.eeg']));
        Fsize=eegfile.bytes;
    otherwise
        error('Bad HDR');
end
Fs=HDR.SampleRate;

%% check if triggers given

%if TT struct not given then obtain it from HDR
if exist('TT','var') == 0  || isempty(TT)
    Trigger= ScouseTom_TrigReadChn(HDR);
    TT= ScouseTom_TrigProcess(Trigger, HDR);
end

%% check expsetup and injection start

if exist('ExpSetup','var') == 0 || isempty(ExpSetup)
    [pathstr,namestr,extstr] = fileparts(HDR.FileName);
    mfilename=fullfile(pathstr,[namestr '_log.mat']);
    %check if _log file is with eeg file
    if exist(mfilename,'file')
        load(mfilename);
    else
        error('Cannot find ExpSetup');
    end
    
end

%find which line in the protocol this starts with
if exist('StartInj','var') ==0 || isempty(StartInj)
    [StartInj] = ScouseTom_data_checkfirstinj(HDR,TT.InjectionSwitches(1,:),ExpSetup.Protocol );
end

Nfreq=size(TT.InjectionSwitches,2);
%find which line in the protocol this starts with
if (exist('FreqNum','var') ==0) || isempty(FreqNum)
    FreqNum=1;   
    if Nfreq > 1
        fprintf(2,'FREQ NOT SELECTED! Choosing 1 for now \n');
    end
end

%% Extract useful stuff
Nelec=size(HDR.InChanSelect,1);
Nprt=size(ExpSetup.Protocol,1);
curWind=TT.InjectionSwitches{FreqNum};

%number of injections found
Ninj=size(curWind,1);

%find number of repeats - rounding up
Nrep=ceil((StartInj(FreqNum)+Ninj)/Nprt);

if Repeat > Nrep   
    error('Repeats too high! Max num is %d\n',Nrep);
end



%% Find relevant injection we wish to view
%get full prt
[prt_full,keep_idx,rem_idx,Elec_inj]=ScouseTom_data_findprt(ExpSetup.Protocol,Nelec);

InjtoView=floor(ChntoView/Nelec)+1; %Injection pair in protocol this corresponds to
CurrentChannel=ChntoView-((InjtoView-1)*Nelec); %which electrode this correspnds to

OtherChannels=1:Nelec;
OtherChannels(CurrentChannel)=[];

InjWind=InjtoView+((Repeat-1)*Nprt); %the data window we want is

InjWind=InjWind+(StartInj(FreqNum)-1); %adjust for starting injection in file

curWind=TT.InjectionSwitches{FreqNum};


%% Find relevant seconds

FirstSample=curWind(InjWind,1);
LastSample=curWind(InjWind,2);

%find the start time in seconds - subtract 1 to give more samples to reduce
%filtering artefacts
StartSec=floor(FirstSample/Fs)-1;
if StartSec < 0
    StartSec =0; %make sure it is >=0
end

%find the stop time in seconds - add 1 to give more samples to reduce
%filtering artefacts
StopSec=ceil(LastSample/Fs)+1;
if StopSec > Nsec
    StopSec=StopSec-1; %remove added second if we go too far
end

Start_Sample=StartSec*Fs;

%% Load and plot data


%plot window is either half an injection either side, or half a second
WindLength=LastSample-FirstSample;

if WindLength > Fs
    FirstPlotSample=floor(FirstSample-Fs/2-Start_Sample);
    LastPlotSample=ceil(LastSample+Fs/2-Start_Sample);
else
    FirstPlotSample=floor(FirstSample-WindLength/2-Start_Sample);
    LastPlotSample=floor(LastSample+WindLength/2-Start_Sample);
end


V=sread(HDR,StopSec-StartSec,StartSec);

%% Plot
t=((0:length(V)-1)./Fs)+StartSec;
injstart=FirstSample./Fs;
injend= LastSample./Fs;


figure
b=colormap(lines);
hold on
%plot other channels

plot(t(FirstPlotSample:LastPlotSample),V(FirstPlotSample:LastPlotSample,OtherChannels),'color',[0.8 0.8 0.8],'Linewidth',0.8)



%plot voltage of interest
plot(t(FirstPlotSample:LastPlotSample),V(FirstPlotSample:LastPlotSample,CurrentChannel),'Linewidth',2,'color',b(1,:))

%plot lines indicating injection switch
plot([injstart,injstart],ylim,'k--','Linewidth',0.8)
plot([injend,injend],ylim,'k--','Linewidth',0.8)


hold off
xlabel('Time (S)');
ylabel('Voltage (uV)');
titlestr=sprintf('Raw Voltage Prt %d/%d Inj %d/%d Elec %d/%d Rep %d/%d',ChntoView,length(prt_full),InjtoView,Nprt,CurrentChannel,Nelec,Repeat,Nrep);

if Nfreq > 1
    titlestr=strcat(titlestr,sprintf(' Freq %d/%d',FreqNum,Nfreq));
end


title(titlestr);



end

