function [] = ScouseTom_TrigView( HDR,Trigger )
%SCOUSETOM_TRIGVIEW Plot trigs in dataset 
%   Detailed explanation goes here


% TO DO:
% needs to not plot repeated values to make it faster

%% get labels from trigger struct if given

trignum=8;

if exist('Trigger') %if trigger file was given, take the labels and the channels which actually have triggers in them

    GoodChn=cellfun(@(x) ~isempty(x),Trigger.RisingEdges);
    TrigDisp=find(GoodChn);
    ChnLabel=(  Trigger.Type(TrigDisp));
    
else %create generic labels if the trigger file is not given
    ChnLabel=({'Chn1','Chn2','Chn3','Chn4','Chn5','Chn6','Chn7','Chn8'});
    TrigDisp=1:trignum;
end


%% Get trigger channel input according to which file type it is
switch HDR.TYPE
    case 'BDF' % biosemi file
        IndicatorPinData=dec2bin(HDR.BDF.Trigger.TYP)-'0';
        IndicatorPinData=IndicatorPinData(:,end-(trignum-1):end); % take only last 8 bits
        IndicatorPinData=fliplr(IndicatorPinData); %sort into LSB
        fname=HDR.FILE.Name;
        TrigPos=HDR.BDF.Trigger.POS;
    case 'EEG'
end

Fs=HDR.SampleRate;

%% Process data

t=TrigPos/Fs;




%% plot that!

figure;

sep=0.5;
hold on
for ichn=(1:length(TrigDisp))

stairs(t,(1.5*(ichn-1)+sep)+IndicatorPinData(:,TrigDisp(ichn)));


end
hold off
ylim([0 1.5*length(TrigDisp)+sep]);
title(['Triggers in dataset: ' fname]);



set(gca,'YTickLabel',ChnLabel,'YTick',1.5*((1:length(TrigDisp))-1)+sep*2,'ygrid','on');
xlabel('Time (s)')
ylabel('Channel');



end

