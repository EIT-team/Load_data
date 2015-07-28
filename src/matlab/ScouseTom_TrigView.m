function [] = ScouseTom_TrigView( HDR,Trigger )
%SCOUSETOM_TRIGVIEW Plot trigs in dataset 
%   Detailed explanation goes here

%% get labels from trigger struct if given

trignum=8;

if exist('Trigger') %if trigger file was given, take the labels and the channels which actually have triggers in them

    GoodChn=cellfun(@(x) ~isempty(x),Trigger.RisingEdges);
    TrigDisp=find(GoodChn);
    ChnLabel=fliplr(  Trigger.Type(TrigDisp));
    
else %create generic labels if the trigger file is not given
    ChnLabel=fliplr({'Chn1','Chn2','Chn3','Chn4','Chn5','Chn6','Chn7','Chn8'});
    TrigDisp=1:trignum;
end


%% Get trigger channel input according to which file type it is
switch HDR.TYPE
    case 'BDF' % biosemi file
        IndicatorPinData=dec2bin(HDR.BDF.ANNONS)-'0';
        IndicatorPinData=IndicatorPinData(:,end-(trignum-1):end); % take only last 8 bits
        IndicatorPinData=fliplr(IndicatorPinData); %sort into LSB
        fname=HDR.FILE.Name;
    case 'EEG'
end

%% plot that!

figure;
strips(IndicatorPinData(:,TrigDisp));
title(['Triggers in dataset: ' fname]);
set(gca,'YTickLabel',ChnLabel);

end

