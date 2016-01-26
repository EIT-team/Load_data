function [] = ScouseTom_TrigView( HDR,Trigger )
%SCOUSETOM_TRIGVIEW Plot trigs in dataset 
%   Detailed explanation goes here



%% Get trigger channel input according to which file type it is
switch HDR.TYPE
    case 'BDF' % biosemi file
        trignum=8;
        [ IndicatorPinData,TrigPos ] = ScouseTom_getbdftrig( HDR,trignum );
        fname=HDR.FILE.Name;
    case 'EEG'
    otherwise
        error('Bad HDR');
end

Fs=HDR.SampleRate;

%% get labels from trigger struct if given

if exist('Trigger') %if trigger file was given, take the labels and the channels which actually have triggers in them

    GoodChn=cellfun(@(x) ~isempty(x),Trigger.RisingEdges);
    TrigDisp=find(GoodChn);
    ChnLabel=(  Trigger.Type(TrigDisp));
    
else %create generic labels if the trigger file is not given
    ChnLabel=({'Chn1','Chn2','Chn3','Chn4','Chn5','Chn6','Chn7','Chn8'});
    TrigDisp=1:trignum;
end


%% Process data

%get time vector
t=TrigPos/Fs;


%% plot that!

figure;
%plot each trigger channel from 0 to 1 with separation of sep. This is like
%plotting using strips, but not shit

sep=0.5;
hold all
for ichn=(1:length(TrigDisp))
stairs(t,(1.5*(ichn-1)+sep)+IndicatorPinData(:,TrigDisp(ichn)));
end
hold off
ylim([0 1.5*length(TrigDisp)+sep]);
title(['Triggers in dataset: ' fname],'interpreter','none');

set(gca,'YTickLabel',ChnLabel,'YTick',1.5*((1:length(TrigDisp))-1)+sep*2,'ygrid','on');
xlabel('Time (s)')
ylabel('Channel');



end

