function [ HDR ] = ScouseTom_getbdfHDR( varargin )
%SCOUSETOM_GETBDFHDR Summary of this function goes here
%   Detailed explanation goes here


%% Ask user for file if not given
%prompt user if no inputs
if isempty(varargin) == 1
    
    [filename, pathname] = uigetfile('*.bdf', 'Choose which Biosemi file with Zcheck to load');
    if isequal(filename,0) || isequal(pathname,0)
        error('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    
    bdfname =fullfile(pathname,filename);

else
    bdfname = varargin{1};
    
end


fprintf('Getting HDR for %s\n',bdfname);

%% check if it exists
if exist(bdfname,'file') ==0
    fprintf(2,'WHOA! FILE DOESNT EXIST!\n');
    HDR=[];
    return
end

%% Read HDR header

% Used in this way we get much less warnings!
%do it once without naming channels to find out how many channels there are
HDR=sopen(bdfname,'r',[],['OVERFLOWDETECTION:OFF','BDF:[4]']);

N_elec=HDR.NS-1; %number of electrodes

%create HDR again, but with correct number of electrodes, Im not exactly
%sure why it doesnt load it properly first
HDR=sopen(bdfname,'r',1:N_elec,['OVERFLOWDETECTION:OFF','BDF:[4]']);


end

