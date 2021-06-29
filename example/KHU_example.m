% Example script to load KHU datasets
% data is stored in the folder selected in a .mat
% first output is corrected BV, second is the KHU structure which contains
% all info about the file, aswell as unscaled data etc.

%% load Baseline

%load boundary voltages of baseline recording, cleaning the data and
%plotting
[BV_baseline, KHUs_baseline]=KHU_Load('../resources/data/KHU/B1',1,1);

%% load Perturbation

%load boundary voltages of baseline recording, cleaning the data but not
%plotting
[BV_perturbation, KHUs_perturbation]=KHU_Load('../resources/data/KHU/P1',1,0);

