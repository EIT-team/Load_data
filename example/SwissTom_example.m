% Example script to load Swisstom dataset
% data is stored in the folder selected in a .mat
% output is structure which contains 
% all info about the file, as well as unscaled data etc.

%% load Baseline

%load boundary voltages of resistor phantom recording. Removing injected
%channels, saturated channels, and plotting data
[STout]=SwissTom_Load('..\resources\data\Swisstom\r-phantom.eit',1,1,1);
