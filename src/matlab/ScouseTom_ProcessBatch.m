function [ output_args ] = ScouseTom_ProcessBatch( dirname )
%SCOUSETOM_PROCESSBATCH Summary of this function goes here
%   Detailed explanation goes here

%% Check or get directory

%user chooses directory where all the .bdfs are
dirname=uigetdir('','Pick the directory where the data is located');
if isempty(dirname)
    error('User Pressed Cancel');
end

%% find all the .bdfs in the directory

files=dir([dirname filesep '*.bdf']);

% check if there are scans actually found
if isempty(files)
    error('No scan files found!');
end

nfiles=length(files);

disp(['Found ' num2str(nfiles) ' .bdf files in directory']);

%ignore small files <1Mb as these are empty 

smallfile=cell2mat({files.bytes})/1e6;
smallfile = smallfile < 1;

if any(smallfile)
    disp(['WARNING! ' num2str(nnz(smallfile)) ' VERY SMALL (<1 Mb) file(s) were detected! These will be ignored']);
    files(smallfile)=[];

end

nfiles=length(files);

%% process each bdf!

tic 
for iFile =1:nfiles
    disp(['Processing file ' num2str(iFile) ' of ' num2str(nfiles) ': ' files(iFile).name]);
    
    ScouseTom_LoadBV(fullfile(dirname,files(iFile).name));
    
    disp('==================================');
    
end

el=toc;

fprintf('All Processing finished in : %.2f seconds\r\n',el);

disp('AWW YISSSS');



end

