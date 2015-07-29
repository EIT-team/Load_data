function [ output_args ] = ScouseTom_LoadZ( HDR,TT,ExpSetup )
%SCOUSETOM_LOADZ Calculates the contact impedance from a two electrode "contact check" 
%   Detailed explanation goes here

%% check inputs are ok here



%% See what type of system we are using 


%the maximum voltage is different for each system, this *should* be in the
%HDR structure somewhere, but I dont know where it is 

switch HDR.TYPE
    case 'BDF' % biosemi file
        
        MaxV=0.5; %500mV range on BioSemi
        

    case 'EEG'
end

%% Get variables

Fs=HDR.SampleRate;

N_elec=ExpSetup.Elec_num; %number of electrodes
N_prt=N_elec; %this is the same as N_elec as the protocol lines is equal to number of electrodes for contact checks



%% Amplitude and freq is FIXED AT THE MOMENT
Amp=1000e-6; 
Frq=1000; 

%% run each impedance check separately

% run each
Z_checks_num=length(TT.InjectionStarts);

for iZ = 1:Z_checks_num


%% read portion of interest 

    %sread needs integer seconds
    Data_start=floor(TT.InjectionStarts(iZ)/Fs);
    Data_end=ceil(TT.InjectionStops(iZ)/Fs);
    Data_length=fix(Data_end-Data_start);
    
    Data_start_s=Data_start*Fs;
    Data_end_s=Data_end*Fs;
    
    %load the data
    disp(['Reading relevant data for contact check :', num2str(iZ)]);
    V=sread(HDR,Data_length,Data_start);
    












end

