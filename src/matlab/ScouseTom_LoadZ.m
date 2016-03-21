function [ Zall ] = ScouseTom_LoadZ( varargin )
%SCOUSETOM_ZLOAD Get Contact Impedance from File selected. This is a
%wrapper for the other TrigProcess and ProcessZ functions
%   Detailed explanation goes here



% Written by Tactile and Tawdry Jimmy 2014, modified 2016


%% Check inputs
%if no inputs, then ask user to pick file with dialogue

if isempty(varargin) == 1
    
    [filename, pathname] = uigetfile({'*.bdf;*.eeg'}, 'Choose which EEG file with Zcheck to load- BioSemi or BrainVision');
    if isequal(filename,0) || isequal(pathname,0)
        error('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    
    fname =fullfile(pathname,filename);
end


if nargin == 1
    fname= varargin{1};
end
if nargin == 2
    
    fname=varargin{1};
    HDR=varargin{2};
    TT=varargin{3};
    
end

if nargin == 3
    fname=varargin{1};
    HDR=varargin{2};
    TT=varargin{3};
    ExpSetup=varargin{4};
    
end

%% Find the file type

[pathstr,namestr,extstr] = fileparts(fname);

%% Load HDR

%if HDR not given then load it
if exist('HDR','var') == 0
    
    HDR=ScouseTom_getHDR(fname);
    
end

%Check HDR
if ~any(strcmp(HDR.TYPE,{'BDF','BrainVision'}))
    error('BAD HDR FILE');
end

%% Get Triggers

%if TT struct not given then obtain it from HDR
if exist('TT','var') == 0
    Trigger= ScouseTom_TrigReadChn(HDR);
    TT= ScouseTom_TrigProcess(Trigger, HDR);
end

%% Get ExpSetup - or find from file

%if ExpSetup not given then load the one. Or if cant be found then just
%take defaults as it doesnt really matter that much!
mfilename=fullfile(pathstr,[namestr '_log.mat']);

if exist('ExpSetup','var') == 0 
    %check if _log file is with eeg file
    if exist(mfilename,'file')
        load(mfilename);
    else
        fprintf(2,'Cannot find ExpSetup, so assuming defaults.this is probably fine\n');
        ExpSetup.Elec_num=HDR.NS-1; %number of electrodes
        ExpSetup.Bad_Elec=[]; %bad electrodes assume none
        ExpSetup.Desc='THIS IS A DUMMY EXPSETUP MADE DURING ZCHECK';
    end
end


%% Process Z

Zall=ScouseTom_ProcessZ(HDR,TT,ExpSetup,1);









end

