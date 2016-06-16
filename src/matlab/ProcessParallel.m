function [dV_signal, t_signal,V_signal,V_baseline,P_signal ,P_baseline,Vmax] = ProcessParallel( fname,BaselineWindow,SignalWindow,TimeStep,InjsSim,BV0,F,plotflag)
%PROCESSPARALLEL Summary of this function goes here
%   Detailed explanation goes here

%%


if exist('plotflag','var') == 0
    plotflag = 1;
end


HDR=ScouseTom_getHDR(fname);

Fs=HDR.SampleRate;
Nchn=size(HDR.InChanSelect,1);

[ StatusChn,TrigPos ] = ScouseTom_geteegtrig(HDR);
TrigPos = TrigPos/Fs ;

Trig_gaps = [ round( diff(TrigPos))];

if isempty(BaselineWindow)
    %if not specified use the first big gap in the triggers
    BaselineWindow(1) = TrigPos(find (Trig_gaps > 1, 1));
    BaselineWindow(2) = TrigPos(find (Trig_gaps > 1, 1)+1);
end

if isempty(SignalWindow)
    %if not specified use the first big gap in the triggers
    SignalWindow(1) = TrigPos(find (Trig_gaps > 1, 1,'last'));
    SignalWindow(2) = TrigPos(find (Trig_gaps > 1, 1,'last')+1);
end



StartSec=max([min([BaselineWindow SignalWindow])-1 0]);
StopSec=max([BaselineWindow SignalWindow])+1;

disp('loading data');
V=sread(HDR,StopSec-StartSec,StartSec);

Vmax= max(V);

if plotflag
    figure
    bar(Vmax)
    xlabel('channel')
    ylabel('uV')
    title('max voltage - USE THIS TO SEE WHAT ELECS YOU WANT TO REMOVE!!!');
    ylim(max(ylim) * [-1 1]);
    drawnow
end

t=(0:length(V)-1)/Fs;

%% Time chunks
baseline_wind=BaselineWindow-StartSec;
signal_wind=SignalWindow-StartSec;


%% Filter Bandwidth

BW=1/TimeStep;

%% Find injection electrodes and freqs

if exist('F','var') == 0
    
    [InjsExp, Freqs] = Find_Injection_Freqs_And_Elecs(V(t<1,:),Fs);
    
    %make sure they are in the same order
    InjsExp=sort(InjsExp,2);
    InjsSim=sort(InjsSim,2);
    
    %put freqs in order to give same protocol as simulation
    [~,Locb]=ismember(InjsSim(:,1),InjsExp(:,1));
    F=Freqs(Locb);
    
end

nFreq=length(F);

%% Demodulate each channel after notch filtering out the other frequencies
Vfull=zeros(size(V,1),size(V,2)*nFreq);
Pfull=zeros(size(Vfull));

for iFreq=1:length(F);
    fprintf('Processing Freq : %d of %d\n',iFreq,length(F));
    cFreq=F(iFreq);
    
    otherfreqs=F;
    otherfreqs(iFreq)=[];
    
    Vf=V;
    
    for inotfreq =1:nFreq-1;
        
        [Bn,An] = butter(2,(otherfreqs(inotfreq)+[-5,5])./(Fs/2),'stop');
        Vf=filtfilt(Bn,An,Vf);
    end
    %notch filter out carriers of other freqs
    
    %band pass filter for demodulation
    [Bdemod,Ademod] = butter(3,(cFreq+[-BW,BW])./(Fs/2));
    [Vdemod,Pdemod]=ScouseTom_data_DemodHilbert(Vf,Bdemod,Ademod);
    
    vidx=(iFreq-1)*Nchn + 1:(iFreq)*Nchn;
    
    Vfull(:,vidx)=Vdemod;
    Pfull(:,vidx)=Pdemod;
end

%put into volts
Vfull=Vfull*1e-6;

%% plot data here
if plotflag
    figure;
    hold on
    h1=plot(StartSec+t,Vfull);
    h2=plot(StartSec+[baseline_wind(1) baseline_wind(1)],ylim,'k--','DisplayName','BaselineWindow');
    h3=plot(StartSec+[baseline_wind(2) baseline_wind(2)],ylim,'k--');
    h4=plot(StartSec+[signal_wind(1) signal_wind(1)],ylim,'r:','DisplayName','SignalWindow');
    h5=plot(StartSec+[signal_wind(2) signal_wind(2)],ylim,'r:');
    hold off
    legend([h2 h4])
    title('Demodulated signals')
    drawnow
end

%% bin data into timesteps to reduce file size

tbins=0:Fs*(TimeStep):length(Vfull)-1;

%put data into bins
[binc,binind]=histc(0:length(Vfull)-1,tbins);
binind(binind ==0)=max(binind+1);

Vbin=zeros(size(tbins,2),size(Vfull,2));
Pbin= zeros(size(Vbin));

for ichn=1:size(Vfull,2)
    tmp=accumarray(binind',Vfull(:,ichn),[],@mean);
    Vbin(:,ichn)=tmp(1:size(tbins,2));
    tmp=accumarray(binind',Pfull(:,ichn),[],@mean);
    Pbin(:,ichn)=tmp(1:size(tbins,2));
end

tbin=(0:length(tbins)-1).*TimeStep;

%change to volts
% Vbin=Vbin*1.e-6;
%take lazy way of finding sign
Vbin=Vbin.*repmat(sign(BV0'),size(Vbin,1),1);

%%
if plotflag
    figure;
    hold on
    h1=plot(StartSec+tbin,Vbin);
    h2=plot(StartSec+[baseline_wind(1) baseline_wind(1)],ylim,'k--','DisplayName','BaselineWindow');
    h3=plot(StartSec+[baseline_wind(2) baseline_wind(2)],ylim,'k--');
    h4=plot(StartSec+[signal_wind(1) signal_wind(1)],ylim,'r:','DisplayName','SignalWindow');
    h5=plot(StartSec+[signal_wind(2) signal_wind(2)],ylim,'r:');
    hold off
    legend([h2 h4])
    title('Binned and signed signal')
    drawnow
end

%% take chunks

V_baseline=mean(Vbin(tbin > baseline_wind(1) & tbin < baseline_wind(2),:));
P_baseline=mean(Pbin(tbin > baseline_wind(1) & tbin < baseline_wind(2),:));

signal_idx=tbin > signal_wind(1) & tbin < signal_wind(2);

V_signal=Vbin(signal_idx,:);
P_signal=Pbin(signal_idx,:);

t_signal=(0:size(signal_idx)-1)/Fs;

%% change in voltage

dV_full=(Vbin-repmat(V_baseline,size(Vbin,1),1));
dV_signal=dV_full(signal_idx,:);

%% filtering comparison

[maxval, maxchn]=max(max(dV_signal));

if plotflag
    figure;
    hold on
    h1=plot(StartSec+tbin,abs(Vbin(:,maxchn)));
    h2=plot(StartSec+t,abs(Vfull(:,maxchn)));
    % h2=plot([baseline_wind(1) baseline_wind(1)],ylim,'k--','DisplayName','BaselineWindow');
    % h3=plot([baseline_wind(2) baseline_wind(2)],ylim,'k--');
    % h4=plot([signal_wind(1) signal_wind(1)],ylim,'r:','DisplayName','SignalWindow');
    % h5=plot([signal_wind(2) signal_wind(2)],ylim,'r:');
    hold off
    % legend([h2 h4])
    title('signal comparison on chn with max dV')
    ylim((maxval*[-1 1])+abs(V_baseline(maxchn)))
    drawnow
end


%% plot final dV
if plotflag
    figure
    hold on
    h1=plot(StartSec+tbin,dV_full);
    ylim(maxval*[-1 1])
    h2=plot(StartSec+[baseline_wind(1) baseline_wind(1)],ylim,'k--','DisplayName','BaselineWindow');
    h3=plot(StartSec+[baseline_wind(2) baseline_wind(2)],ylim,'k--');
    h4=plot(StartSec+[signal_wind(1) signal_wind(1)],ylim,'r:','DisplayName','SignalWindow');
    h5=plot(StartSec+[signal_wind(2) signal_wind(2)],ylim,'r:');
    hold off
    legend([h2 h4])
    title('Voltage Change whole data set')
    drawnow
end

end

