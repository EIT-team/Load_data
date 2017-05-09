% Example script to load Swisstom dataset
% data is stored in the folder selected in a .mat
% output is structure which contains 
% all info about the file, as well as unscaled data etc.

%% load Baseline

%load boundary voltages of resistor phantom recording. Removing injected
%channels, saturated channels, and plotting data
[BVstruc]=ScouseTom_LoadBV('..\resources\data\ScouseTom\forreal_2_clean.bdf');

%% load Baseline

%load boundary voltages of resistor phantom recording. Removing injected
%channels, saturated channels, and plotting data
[BVstruc]=ScouseTom_LoadBV('..\resources\data\ScouseTom\MultiFreq.bdf');