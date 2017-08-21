function [ HDR ] = ScouseTom_getHDR( varargin )
%SCOUSETOM_GETBDFHDR Summary of this function goes here
%   Detailed explanation goes here

MinFileSize=1e6; %Minimum file size in bytes

%% Ask user for file if not given
%prompt user if no inputs

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

fprintf('Getting HDR for %s\n',fname);

%% check if it exists
if exist(fname,'file') ==0
    fprintf(2,'WHOA! FILE DOESNT EXIST!\n');
    HDR.TYPE='NULL';
    return
end


%% Find the file type

[pathstr,namestr,extstr] = fileparts(fname);

%%
%use function for correct file type
switch extstr
    case '.bdf'
        % Used in this way we get much less warnings!
        %do it once without naming channels to find out how many channels there are
        HDR=sopen(fname,'r',[],['OVERFLOWDETECTION:OFF','BDF:[4]']);
        N_elec=size(HDR.InChanSelect,1); %number of electrodes
        %create HDR again, but with correct number of electrodes, Im not exactly
        %sure why it doesnt load it properly first
        HDR=sopen(fname,'r',1:N_elec,['OVERFLOWDETECTION:OFF','BDF:[4]']);
    case {'.eeg','.vhdr','.vmrk'}
        %This is much simpler
        HDR=sopen(fname,'r',[],['OVERFLOWDETECTION:OFF']);
    otherwise
        HDR=[];
        error('Unknown file type');
end


%% Check HDR

%Check HDR is in either of expected formats
if ~any(strcmp(HDR.TYPE,{'BDF','BrainVision'}))
    HDR=[];
    error('BAD HDR FILE');
end


switch HDR.TYPE
    case 'BDF' % biosemi file
        %biosemi is stored as 1sec long records. so number of seconds is
        %Nrec
        Nsec=HDR.NRec;
        %biosemi file size is stored in bytes
        Fsize=HDR.FILE.size;
    case 'BrainVision'
        Nsec=ceil(HDR.SPR/HDR.SampleRate);
        %actichamp HDR points to header .vhdr file, so filesize is wrong
        eegfile=dir(fullfile(HDR.FILE.Path,[HDR.FILE.Name '.eeg']));
        Fsize=eegfile.bytes;
    otherwise
        error('Bad HDR');
end

if Fsize < MinFileSize
    fprintf(2,'WHOA! FILE IS WAY TOO SMALL! DID YOU FORGET TO START RECORDING?\n');
    HDR=[];
    return
end



end

