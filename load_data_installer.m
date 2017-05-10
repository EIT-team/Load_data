%% Install Biosig Module

homedir=pwd;

cd('./lib/biosig4octmat-3.0/')
biosig_installer

cd(homedir);

%% Add matlab dirs to path

addpath(genpath([pwd filesep 'src' filesep 'matlab']));

%% 

%savepath

disp('Load data path set ok');