% installer for Load_data repository. Loads this modified version of Biosig library then adds the src folder
%% Install Biosig Module

homedir=pwd;

cd('./lib/biosig4octmat-3.0/')
biosig_installer

cd(homedir);

%% Add matlab dirs to path

addpath(genpath([pwd filesep 'src' filesep 'matlab']));

%%

%savepath
disp('###################################');
disp('Load data path set ok. Run SAVEPATH to save these changes');
disp('Please run the examples in the ./examples folder to check everything is working')
