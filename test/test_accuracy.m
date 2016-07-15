% script to test the accuracy of the processing code

Fs=16384;
BW=50;
chn=5;

fc = 500;

injtime=2;
num_inj = 2;

inj = [ 2 4];
Amp_Inj = 500;
Amp_Meas = 150;


s_cyc = Fs/fc;

MeasPhase = 45/360;

MeasDelay = ceil(MeasPhase * s_cyc);

MeasPhaseActual= (MeasDelay/s_cyc);



Totaltime=injtime*num_inj;
t=(0:Totaltime*Fs)/Fs;
%make sin wave
v= sin(2*pi*fc*t);
% change amplitude
v_m=v *Amp_Meas;
v_m = circshift(v_m,[1,MeasDelay]);

% change amplitude
v_i=v *Amp_Inj;


%%
V=repmat(v_m,chn,1);

V(inj(1),:) = v_i;
V(inj(2),:) = v_i;

V=V';

InjectionWindows =[1 (length(V)-10)];



%%
%find the corresponding filter settings
[trim_demod,B,A,Fc]=ScouseTom_data_GetFilterTrim(V(:,inj(1)),Fs,BW,0 );

[ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B,A);


%%
[Vmag,PhaseRaw,VmagSTD,PhaseRawSTD]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,trim_demod,InjectionWindows);

[Phase]=ScouseTom_data_PhaseEst(PhaseRaw,inj,1);

