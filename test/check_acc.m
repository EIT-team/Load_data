function [ Amp_error, Phase_error] = check_acc( Fc,InjTime,Amp_Inj,Amp_Meas,InjPhase,MeasPhaseDiff )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%
num_inj = 1;

inj = [ 2 4];

Fs=16384;
BW=50;
chn=5;

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


Totaltime=InjTime*num_inj;
t=(0:Totaltime*Fs)/Fs;
%make sin wave

v_m= Amp_Meas*sin(2*pi*Fc*t+(pi*MeasPhase/180));

% change amplitude
v_i=Amp_Inj*sin(2*pi*Fc*t+(pi*InjPhase/180));


%%
V=repmat(v_m,chn,1);

V(inj(1),:) = v_i;
V(inj(2),:) = v_i;

V=V';

InjectionWindows =[1 (length(V)-10)];

%%
%find the corresponding filter settings
[trim_demod,B,A,Fc_found]=ScouseTom_data_GetFilterTrim(V(:,inj(1)),Fs,BW,0 );

[ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B,A);


%%
[Vmag,PhaseRaw,VmagSTD,PhaseRawSTD]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,trim_demod,InjectionWindows);

[Phase]=ScouseTom_data_PhaseEst(PhaseRaw,inj,1);

Phase_deg=(180/pi) * (Phase);

Phase_error = PhaseActual - Phase_deg';

Amp_error = AmpActual - Vmag';

fprintf('Amp error : %.6f, Phase error : %.6f\n',mean(Amp_error),mean(Phase_error));
end

