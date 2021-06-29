
Fc=20;

Fcur = 15;
FreqNum = size(Fcur,2);


NumInj= 5;

Cycles = 128;
T=(1./Fcur); %Period in s
InjTime=(T.*Cycles)*NumInj;
InjSamples=(T.*Cycles)*Fs;

% InjTime=10;

Amp_Inj = 500;
Amp_Meas = 150;
InjPhase=0;
MeasPhaseDiff=-30;

inj = [ 2 4];

DCoffset=10+10*(1:NumInj);
DCoffsetInj=100 + 50*(1:NumInj);


Fs=16384;
BW=100;
chn=5;

Padsec=4;

%% Create ideal values, and voltages

AmpActual= repmat(Amp_Meas,chn,1);
AmpActual(inj)=Amp_Inj;

a=MeasPhaseDiff;

%force normal valus of deg
a = mod(a,360); % [0 2pi)

% shift
jj = a > 180;
a(jj) = a(jj) - 360;
jj = a < 0 - 180;
a(jj) = a(jj) + 360;

MeasPhaseDiff_corr = a;

MeasPhase = InjPhase + MeasPhaseDiff_corr; %deg
PhaseActual= repmat(MeasPhaseDiff_corr,chn,1);
PhaseActual(inj)=0;

Totaltime=ceil(InjTime + 2*Padsec);
t = 0:1/Fs:Totaltime-1/Fs;
%make sin wave

v_m= Amp_Meas*sin(2*pi*Fc*t+(pi*MeasPhase/180));

% change amplitude
v_i=Amp_Inj*sin(2*pi*Fc*t+(pi*InjPhase/180));

%%
V=repmat(v_m,chn,1);

V(inj(1),:) = v_i;
V(inj(2),:) = v_i;

V=V';

%pad with a second of data either side, so the hilbert is more realistic
datastart = round(Padsec*Fs);
dataend = round((Padsec+InjTime)*Fs);

V(1:datastart,:)=V(1:datastart,:)*0.1;
V(dataend:end,:)=V(dataend:end,:)*0.1;

InjectionWindows(:,1) = (datastart + ceil(InjSamples * (0:NumInj-1)))';
InjectionWindows(:,2) = [InjectionWindows(2:end,1) ; dataend];


DCoffsetVec = ones(chn,1)*DCoffset;
DCoffsetVec(inj(1),:) = DCoffsetInj;
DCoffsetVec(inj(2),:) = DCoffsetInj;

for iInj = 1:NumInj

V(InjectionWindows(iInj,1):InjectionWindows(iInj,2),:)=bsxfun(@plus,V(InjectionWindows(iInj,1):InjectionWindows(iInj,2),:),DCoffsetVec(:,iInj)');

end


%% 

Voff =mean(V(datastart:dataend,:));

plot(bsxfun(@minus,V,Voff))







