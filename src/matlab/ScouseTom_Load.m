function [ OutStruc ] = ScouseTom_Load( varargin )
% [OutStruc] = ScouseTom_Load( varargin )
%ScouseTom_LoadBV Demodulate the boundary voltages in a EEG file created by
%a ScouseTom system. This is essentially a wrapper for the other
%TrigReadChn and TrigProcess and ProcessBV or ProcessZ functions.
%
% First it gets the HDR structure of the .bdf or .eeg file -
% prompting to choose is not specified by fname or HDR.
%
% Second the Trigger channels are processed to find if this is a "normal"
% injection or a contact check - USer can give TT structure directly
%
% Demodulation is then performed on this dataset either by calling
% ScouseTom_ProcessBV for conventional EIT recordings, or
% ScouseTom_ProcessZ for contact checks, with the results graphed in this
% case.
%
% Example uses
% [] = ScouseTom_LoadBV() - prompts user to select file through GUI
%
% [] = ScouseTom_LoadBV(fname) - loads file given by fname e.g.
%
% ScouseTom_LoadBV('../../resources/data/ScouseTom/forreal_2_clean.bdf');
%
% [] = ScouseTom_LoadBV(HDR) - loads file for given HDR, useful to avoid
% having to calculate the HDR repeatedly when changing demod settings
%
% HDR = ScouseTom_getHDR('../../resources/data/ScouseTom/forreal_2_clean.bdf');
% ScouseTom_LoadBV(HDR);
%
% Written by the avuncular yet bohemian Jimmy 2016

%% Check inputs

if nargin >= 1
    if ischar(varargin{1})
        
        fname= varargin{1};
        
    elseif isstruct(varargin{1})
        HDR = varargin{1};
    end
else
    %if no inputs, then ask user to pick file with dialogue
    [filename, pathname] = uigetfile({'*.bdf;*.eeg'}, 'Choose which EEG file to load- BioSemi or BrainVision');
    if isequal(filename,0) || isequal(pathname,0)
        error('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    fname =fullfile(pathname,filename);
end

if nargin >= 2
    TT=varargin{2};
end

if nargin >= 3
    ExpSetup=varargin{3};
end

if nargin >= 4
    PlotFlag=varargin{4};
end


%% Load HDR

%if HDR not given then load it
if exist('HDR','var') == 0 || isempty(HDR)
    HDR=ScouseTom_getHDR(fname);
end

%Check HDR is legit
if ~any(strcmp(HDR.TYPE,{'BDF','BrainVision'}))
    error('BAD HDR FILE');
end

%% Get Triggers

%if TT struct not given then obtain it from HDR
if exist('TT','var') == 0 || isempty(TT)
    Trigger= ScouseTom_TrigReadChn(HDR);
    TT= ScouseTom_TrigProcess(Trigger, HDR);
end

%% Find log file

[pathstr,namestr,extstr] = fileparts(HDR.FileName);
%if ExpSetup not given then load the one matching the filename (saved
%by ScouseTom_Start
mfilename=fullfile(pathstr,[namestr '_log.mat']);

%% Process normal injection stuff
% Starts of normal protocols should be saved in TT structure
if ~isempty(TT.InjectionStarts)
    % Get ExpSetup - or find from file
    if exist('ExpSetup','var') == 0 || isempty(ExpSetup)
        %check if _log file is with eeg file
        if exist(mfilename,'file')
            load(mfilename);
        else
            error('Cannot find ExpSetup');
        end
    end
    %process the boundary voltages and only take the output structure
    [OutStruc]=ScouseTom_ProcessBV(HDR,TT,ExpSetup);
end

%% Process contact impedance checks - if any found

if ~isempty(TT.Contact.InjectionStarts)
    
    fprintf(2,'File contains Z check. Calling ScouseTom_ProcessZ\n');
    
    if exist('ExpSetup','var') == 0 || isempty(ExpSetup)
        %check if _log file is with eeg file
        if exist(mfilename,'file')
            load(mfilename);
        else % use defaults if not-this only pertains to skipped channels during contact check and number of electrodes
            fprintf(2,'Cannot find ExpSetup, so assuming defaults\n');
            %             ExpSetup.Elec_num=HDR.NS-1; %number of electrodes
            %             ExpSetup.Bad_Elec=[]; %bad electrodes assume none
            %             ExpSetup.Desc='THIS IS A DUMMY EXPSETUP MADE DURING ZCHECK';
        end
    end
    
    if exist('PlotFlag','var') ==0 || isempty(PlotFlag)
        PlotFlag =1;
    end
    
    
    % Process Z check, outputing figures
    OutStruc=ScouseTom_ProcessZ(HDR,TT,[],PlotFlag);
    
end
end

