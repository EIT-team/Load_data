function [ STout,SwissTom_raw] = SwissTom_Load( fname,rem_flag,sat_flag,plot_flag )
% [ STout,SwissTom_raw] = SwissTom_Load( fname,rem_flag,sat_flag,plot_flag )
%   loads data from swisstom .eit file into UCL/view_data format. Extra output
%   is the raw output of the swisstom loader (you know, for safe keeping)
%   loads eit data set from swisstom system, saves in .mat
%   file the voltages across frames format. Also makes the protocol for the given data
%   set
%
%   Data is loaded and a protocol file created based on the offset 
%   (i.e. the spacing between injection and measurement electrodes) value.
%   User is then prompted if they want to remove the measurements from
%   injection electrodes. A new protocol file is saved if this is the case

%   inputs - filename - uses the same filename to load and save .mat and
%   .prt file
%   Clean flag - 1 to remove injection channels otherwise 0
%   Plot flag - 1 to plot data after loading otherwise 0
%
%   outputs: -
%   STout - the converted values in V R X Z Theta
%   SwissTom_raw - raw output from swisstom loader
%
%   files stored:
%   fname-protocol.prt - protocol file generated from settings in eit file
%   fname-SWISSTOMDATA - Structure with all the data in it
%

%% load file

%prompt user if no inputs
if exist('fname','var') == 0
    load_flag=1;
elseif isempty(fname) ==1
    load_flag=1;
else
    load_flag=0;
end


if load_flag ==1
    % ask user for file path if we dont have one given
    [filename, pathname] = uigetfile('*.eit', 'Choose which SwissTom file to load');
    if isequal(filename,0) || isequal(pathname,0)
        error('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    
    fname =fullfile(pathname,filename);
end

disp('Loading SwissTom file');
%read data using swisstom matlab function
SwissTom_raw=SwissTom_EIT_reader(fname);

%% file info
N_elec=SwissTom_raw.tic.num_elec;
offsetvalue=SwissTom_raw.tic.injectionPattern;
N_frames=SwissTom_raw.numberOfFrames;
ScanTime=(SwissTom_raw.finalTimestamp-SwissTom_raw.initialTimestamp)/1000; %time taken in seconds

disp('---------------------------------');
fprintf('Frequency Detected: %d kHz\n',SwissTom_raw.tic.currentFreq);
fprintf('Current Injected: %d mA RMS OR PEAK TO PEAK THERE IS CONFUSION HERE \r',SwissTom_raw.tic.injectionCurrent);
fprintf('Number of frames/repeats/scans : %d which took %.2f seconds\r',N_frames,ScanTime);
fprintf('Number of electrodes: %d IF THIS ISNT 32 THEN SOMETHINGS WRONG \r',N_elec);
disp('---------------------------------');

disp('Now lets crack on...');

%get just the filename and path from this
[pathstr, namestr]=fileparts(fname);

%% voltage measurments - "frames"

disp('Getting data from frames and converting into mV');

%get the data from the .eit file and convert to mV and Ohms
STout=IQtoVZ(SwissTom_raw);


%% remove injecting electrodes if user wants

if exist('rem_flag','var') ==0
    
    rem_ans=questdlg('Remove the measurements on the injecting electrodes? The answer is probably yes','Ignore injecting electrode measurements','Yes, that seems entirely sensible','No - caution to the wind!','Yes, that seems entirely sensible');
    
    if strcmp(rem_ans,'Yes, that seems entirely sensible') == 1
        disp('Getting rid of the voltages recorded on injection channels');
        rem_flag=1;
        
        disp('Writing new protocol file with injection measurements removed');
        
    else
        disp('User chose not to remove injection electrodes');
        rem_flag=0;
    end
    
end
%% protocol file

fprintf('Injection pattern with offset %d detected. Making .prt file using this value...',offsetvalue);

[SwissTomProtocol, rem_idx,keep_idx]=SwissTom_MakeProtocol(offsetvalue,rem_flag,fullfile(pathstr,[namestr,'-protocol.prt']));

fprintf('done!\n');


%% remove saturated channels?

satfound=nnz(STout.sat);

if exist('sat_flag','var') ==0
    
    if satfound >0
        
        sat_ans=questdlg(['Warning! ' num2str(satfound) ' saturated channels found. Add these to the rem_idx and remove from keep_idx?'],'Remove Saturated channels?','Yes Please!','No, I like bad data','Yes Please!');
        
        if strcmp(sat_ans,'Yes Please!') ==1
            sat_flag=1;
            disp('User chose to remove saturated channels');
            
        else
            disp('User chose NOT to remove saturated channels');
            sat_flag=0;
        end
    end
    
end

if sat_flag==1
    
    %find the new channels to remove
    satchn=find(STout.sat);
    %find unique values as most likely duplicates
    rem_idx=unique([rem_idx satchn]);
    keep_idx=1:length(STout.sat);
    keep_idx(rem_idx)=[];
    
end
%% store in strucutre

disp('Sorting out which voltages are needed');

%store data in raw and clean structures
STout.rem_idx=rem_idx;
STout.keep_idx=keep_idx;
STout.filename=fname;

STout.raw.prt=SwissTomProtocol;
STout.freq=(SwissTom_raw.tic.currentFreq)*1000;

%store the BVs
STout.BV_full=(STout.raw.v);
STout.BV=(STout.raw.v(keep_idx,:));
STout.prt_full=SwissTomProtocol;
STout.prt=SwissTomProtocol(keep_idx,:);
disp('All done!');


%% save dat shit
save(fullfile(pathstr,[namestr,'-SWISSTOMDATA.mat']),'-struct','STout');

%% Plot data

if exist('plot_flag','var') ==0
    
    plot_flag=questdlg('Plot data?','Do you wanna plot the data?','Yes Please','Nah mate','Yes Please');
    
    if strcmp(plot_flag,'Yes Please') == 1
        plot_flag =1;
    else
        plot_flag =0;
    end
    
end


if plot_flag == 1
    
    disp('Rendering data and outputting to pixel buffer for human fuzzy logic occular assessment');
    Voltages=STout.BV;
    
    figure;
    plot(mean(abs(Voltages),2))
    xlabel('Electrode Combination');
    ylabel('|Voltage| mV');
    title(['Averaged |Voltage| in file: ', namestr]);
    
    figure
    plot(real(Voltages'))
    xlabel('Frame/Scan/Repeat');
    ylabel('|Voltage| V');
    title(['Voltage across measurements in file: ', namestr]);
    
    
    noiseest=abs((100*( std(Voltages,0,2))./(mean(Voltages,2))));
    figure
    hist(noiseest,100)
    ylabel('Number of channels');
    xlabel('% Noise - Std Dev / Mean');
    title(['Noise estimate in in file: ', namestr]);
    
    figure;
    plot(abs(mean(Voltages,2)),abs(std(Voltages,0,2)),'x')
    xlabel('Mean Voltage mV')
    ylabel('Abs Standard Deviation')
    title(['Noise against boundary voltage in file: ', namestr]);
    
    figure;
    plot(abs(mean(Voltages,2)),100*(abs(std(Voltages,0,2))./abs(mean(Voltages,2))),'x')
    xlabel('Mean Voltage mV')
    ylabel('StdDev/Mean %of mean')
    %     set(gca,'YScale','log')
    title(['Normalised Noise against boundary voltage in file: ', namestr]);
    
end

end


function [out]=IQtoVZ(SwissTom_raw)
% this is the new conversion to voltage and ohms based on the pseudocode
% from peter krammer. 

%pga gains
pga0gains = [1 2 5 10];
pga1gains = [1, 10, 20, 30, 40, 60, 80, 120, 157, 0.2];

%factors
ad_factor = 2^15;
dac1NumberOfSamples = SwissTom_raw.sbc.nsmp;
dac0NumberOfSamples = SwissTom_raw.sbc.nsmp;

%get values from EIT file
pga0gain=pga0gains(SwissTom_raw.sbc.pga0_g+1);
pga1gain=pga1gains(SwissTom_raw.sbc.pga1_g+1);
current=SwissTom_raw.tic.injectionCurrent/1000; %current given in mA

%calculate gain factors - this is where the unknown factors come in
diff_gain= 2 * 0.939167 * pga0gain *pga1gain;
abs_gain = 0.2074 * 0.83410716;

%factor to scale IQ values by
iq_factor = dac0NumberOfSamples /  4^(floor(log(dac0NumberOfSamples)/log(4)));

%conversion factor to scale CORRECT iq values by to get volts and ohms
amplitudefactor= 0.00025*4 / diff_gain;
impedancefactor= (0.00025*4) /(abs_gain * current);
%used to get phase
phasefactor=180/pi;

%load data from frame at a time - IQ values are stored as the digital
%values from demodulation. So first convert them into a voltages
for ii = 1:length(SwissTom_raw.frames)
    %voltage data
    IQ.I(ii,:) = SwissTom_raw.frames(ii).iqPayload(1:2:end)/iq_factor / ad_factor;
    IQ.Q(ii,:) = SwissTom_raw.frames(ii).iqPayload(2:2:end)/iq_factor / ad_factor;
    %Impedance data
    IQ.IZ(ii,:) = SwissTom_raw.frames(ii).voltageInjectionPayload(1:2:end)/iq_factor / ad_factor;
    IQ.QZ(ii,:) = SwissTom_raw.frames(ii).voltageInjectionPayload(2:2:end)/iq_factor / ad_factor;
end

%get the magnitude
v =  sqrt(IQ.I.^2 + IQ.Q.^2)*amplitudefactor;

%check for saturation
sat_thres=4095; %value read off graph in swisstom software, this needs confirming
IQmag=sqrt(IQ.I.^2 + IQ.Q.^2);
IQthres=IQmag > sat_thres;

sat=any(IQthres);

if any(sat,2)
    disp(['Warning! ' num2str(nnz(sat)) ' saturated channels found :(']);
end

%get the phase of the voltage
theta=atan(IQ.Q./IQ.I)*phasefactor;

%convert into real and imaginary components
R=v.*cos(theta);
X=v.*sin(theta);

%calculate contact impedance
z = sqrt(IQ.IZ.^2 + IQ.QZ.^2)*impedancefactor;

out.raw.v=v';
out.raw.z=z';
out.raw.R=R';
out.raw.X=X';
out.raw.theta=theta';
out.raw.IQthres=IQthres;
out.sat=sat;
out.raw.factors.impedacefactor=impedancefactor;
out.raw.factors.phasefactor=phasefactor;
out.raw.factors.IQtoRawV=ad_factor;
out.raw.factors.RawVtoActualV=amplitudefactor;

end


function [ SwissTomProtocolFull, rem_idx,keep_idx ] = SwissTom_MakeProtocol( offsetvalue,removeflag,varargin )
%   SwissTom_MakeProtocol makes .prt file for given swisstom "offset" value
%
%   Swisstom system has a fixed measurement and injection offset.
%   This means the number of electrodes between the source and sink
%   - so an offset of 0 will inject between 1 and 2 and also measure
%   between 1 and 2. an offset of 3 would injetion between 1 and 5 and also
%   measure between 1 and 5. 
%
%   inputs:
%   offsetvalue - setting in the swisstom software. determines which
%   electrodes are used in which lines of the protocol
%   [filename]= filename to save .prt file
%   outputs:
%   SwissTomProtocol - matrix of protocol

%% check inputs

if ceil(offsetvalue) ~= floor(offsetvalue)
    warning('non integer offset, rounding down');
    offsetval = floor(offsetvalue);
end

if offsetvalue < 0
    error('offset must be greater or equal to 0');
elseif offsetvalue >14
    error('offset must be less than 14');
end

%% file dialog

if isempty(varargin) == 1
    
    [filename, pathname] = uiputfile('*.prt', 'Choose where to save the .prt file');
    if isequal(filename,0) || isequal(pathname,0)
        error('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    
    fname =fullfile(pathname,filename);
else
    fname=varargin{1};
end

%% make dat protocol motherfucker!!!

%injection pairs
inj_p=1:32; %sources
inj_n=circshift(inj_p,[1 -(offsetvalue+1)]); %sinks (shifted by offset+1)

%measurementpairs - kept separate because reasons
meas_p=1:32;%sources
meas_n=circshift(meas_p,[1 -(offsetvalue+1)]);%sinks (shifted by offset+1)

N_prt=length(inj_p); %yes these will never change and always be 32, but hey
N_elec=length(meas_p);


%populate matrix of protocol lines
for iPrt = 1:N_prt;
    
    idx=((iPrt-1)*N_elec)+1;
    
    SwissTomProtocol(idx:idx+N_elec-1,:)=[repmat([inj_p(iPrt) inj_n(iPrt)],N_elec,1) meas_p' meas_n'];
end
%save the entire protocol;
SwissTomProtocolFull=SwissTomProtocol;

%% remove data?

if removeflag == 1
    rem_idx=[];
    for iPrt = 1:size(SwissTomProtocol,1)
        if any(ismember(SwissTomProtocol(iPrt,1:2),SwissTomProtocol(iPrt,3:4))) ==1
            
            rem_idx=[rem_idx,iPrt];
        end
    end
    keep_idx=setdiff(1:length(SwissTomProtocol),rem_idx);
    SwissTomProtocol(rem_idx,:)=[];
    
    s='# SwissTom specific protocol with OffsetValue %d. Measurements on injection electrodes removed. Automagically created by the winsome jimmy\r\n';
else
    rem_idx=[];
    keep_idx=1:length(SwissTomProtocol);
    s='# SwissTom specific protocol with OffsetValue %d. Automagically created by the winsome jimmy\r\n';
    
end

%% write data protocol!!!

%write them to a file with header
fid=fopen(fname,'w+');

fprintf(fid,s,offsetvalue);

for iThing =1:length(SwissTomProtocol)
    fprintf(fid,'%d %d %d %d \r\n',SwissTomProtocol(iThing,1),SwissTomProtocol(iThing,2),SwissTomProtocol(iThing,3),SwissTomProtocol(iThing,4));
end

fclose(fid);


end


