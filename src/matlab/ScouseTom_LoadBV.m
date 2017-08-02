function [ BVstruc ] = ScouseTom_LoadBV( varargin )
%SCOUSETOM_ZLOAD Get Contact Impedance from File selected. This is a
%wrapper for the other TrigProcess and ProcessZ functions
%   Detailed explanation goes here


% Written by the avuncular yet bohemian Jimmy 2016

%% Check inputs
%if no inputs, then ask user to pick file with dialogue

if isempty(varargin) == 1
    
    [filename, pathname] = uigetfile({'*.bdf;*.eeg'}, 'Choose which EEG file to load- BioSemi or BrainVision');
    if isequal(filename,0) || isequal(pathname,0)
        error('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    
    fname =fullfile(pathname,filename);
    
else
    fname = varargin{1};
    
end

if nargin >= 1
    fname= varargin{1};
end
if nargin >= 2
    HDR=varargin{2};
end

if nargin >= 3
    TT=varargin{3};
end

if nargin >= 4
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

%% Do Normal injection stuff

if ~isempty(TT.InjectionStarts)
    %% Get ExpSetup - or find from file
    
    %if ExpSetup not given then load the one
    
    mfilename=fullfile(pathstr,[namestr '_log.mat']);
    
    if exist('ExpSetup','var') == 0
        
        %check if _log file is with eeg file
        if exist(mfilename,'file')
            load(mfilename);
        else
            error('Cannot find ExpSetup');
        end
        
    end
    
    
    %% Process Boundary Voltages
    
    %process the boundary voltages and only take the output structure
    [BVstruc]=ScouseTom_ProcessBV(HDR,TT,ExpSetup);
    
    
    
    
    
end



if ~isempty(TT.Contact.InjectionStarts)
    fprintf(2,'File contains Z check. Calling ScouseTom_LoadZ\n');
    
    [ BVstruc ] = ScouseTom_LoadZ(fname, HDR,TT );


end



end

