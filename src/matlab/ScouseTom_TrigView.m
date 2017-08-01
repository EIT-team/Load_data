function [] = ScouseTom_TrigView( HDR,XaxisSampleFlag,Trigger)
% [] = ScouseTom_TrigView( HDR,XaxisSampleFlag,Trigger )
% SCOUSETOM_TRIGVIEW Plots the triggers across all digital channels in
% dataset, in something resembling strips. Useful in debugging processing
% [] = ScouseTom_TrigView( HDR,XaxisSampleFlag,Trigger )
% Inputs:
% HDR - HDR from ScouseTom_getHDR
% XaxisSampleFlag (optional) - Choose samples or seconds for X axis : 0 = sec [def], 1 = samples
% Trigger (optional) - Provide output of ScouseTom_TrigReadChn to label
% channels 


%% Get trigger channel input according to which file type it is
switch HDR.TYPE
    case 'BDF' % biosemi file
        trignum=8;
        [ StatusChns,TrigPos ] = ScouseTom_getbdftrig( HDR,trignum );
        fname=HDR.FILE.Name;
        plotstr='-';
    case 'BrainVision'
        [ StatusChns,TrigPos ] = ScouseTom_geteegtrig( HDR );
        fname=HDR.FILE.Name;
        trignum=size(StatusChns,2);
        plotstr='-*';
    otherwise
        error('Bad HDR');
end

Fs=HDR.SampleRate;

%% get labels from trigger struct if given

if exist('Trigger','var') %if trigger file was given, take the labels and the channels which actually have triggers in them
    
    GoodChn=cellfun(@(x) ~isempty(x),Trigger.RisingEdges);
    TrigDisp=find(GoodChn);
    ChnLabel=(  Trigger.Type(TrigDisp));
    
else %create generic labels if the trigger file is not given
    ChnLabel=({'Chn1','Chn2','Chn3','Chn4','Chn5','Chn6','Chn7','Chn8'});
    TrigDisp=1:trignum;
end

if ~exist('XaxisSampleFlag','var')
    XaxisSampleFlag=0;
end

%% Process data

%get time vector
t=TrigPos/Fs;


%% plot that!

figure;
if XaxisSampleFlag
    xdata=TrigPos;
    xlabel('Sample')
    
else
    xdata=t;
    xlabel('Time (s)')
end

%plot each trigger channel from 0 to 1 with separation of sep. This is like
%plotting using strips, but actually works 

sep=0.5;
hold all
for ichn=(1:length(TrigDisp))
    stairs(xdata,(1.5*(ichn-1)+sep)+StatusChns(:,TrigDisp(ichn)),plotstr);
end
hold off
ylim([0 1.5*length(TrigDisp)+sep]);
title(['Triggers in dataset: ' fname],'interpreter','none');

set(gca,'YTickLabel',ChnLabel,'YTick',1.5*((1:length(TrigDisp))-1)+sep*2,'ygrid','on');

ylabel('Channel');



end

