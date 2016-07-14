% script to test the accuracy of the processing code

Fs=16384;
Bandwidth=50;
chn=10;

fc = 1000;

injtime=2;
num_inj = 2;

Totaltime=injtime*num_inj;
t=(0:Totaltime*Fs)/Fs;

v= sin(2*pi*fc*t);

Amp_Inj = 500;
Amp_Meas = 150;

Amp





% V=




%%
%find the corresponding filter settings
[trim_demod,B,A,Fc]=ScouseTom_data_GetFilterTrim(V(tmpidx),Fs,BandWidth,0 );

[ Vdata_demod,Pdata_demod ] = ScouseTom_data_DemodHilbert( V,B,A);

[Vmag{iFreq}(:,curChn),PhaseRaw{iFreq}(:,curChn),VmagSTD{iFreq}(:,curChn),PhaseRawSTD{iFreq}(:,curChn)]=ScouseTom_data_getBV(Vdata_demod,Pdata_demod,trim_demod,InjectionWindows{iFreq}-Start_Sample);


