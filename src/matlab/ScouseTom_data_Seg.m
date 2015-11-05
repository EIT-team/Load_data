function [ Vseg,Pseg, lastPrt ] = ScouseTom_data_Seg( Vdata,Pdata,InjSwitches,FreqStartSwitches,FreqStopSwitches,t_rem,N_prt,N_elec,Fs,varargin )
%segment_ind segments data stream based on timepoints InjSwitches - taken
%from HDR. data is removed hereto remove switching artefact, but this is
%not used that much. most data removal is done in demodulation
%   Detailed explanation goes here

%% data info

if isempty(varargin) ==1
    iPrt=1;
else
    tmp=varargin{1};
    %make sure it is a number
    if isnumeric(tmp) ==1
        %this is unessesary now
        if tmp ~= 1
            %if the next line is bigger than the number of protocol lines
            %then "wrap" it round
            if tmp>N_prt
                iPrt=tmp-N_prt;
            else
                iPrt=tmp;
            end
        else
            iPrt=1;
        end
    end
end


%% info for segmenting

s_rem=fix(t_rem*Fs); %samples to remove

N_switch=size(InjSwitches,1); %number of switches

N_rep=ceil(N_switch/N_prt); %number of repeats of protocol

%are we taking a specific bit within this injection for multiFreqMode?
if (~isempty(FreqStartSwitches) && ~isempty(FreqStopSwitches))
    MultiFreqMode =1;
    %then the interval is the biggest gap between freqstart and freqstop
    s_int=max(FreqStopSwitches-FreqStartSwitches);
    
    %calculate vector of trimmed segments - reference everything to the
    %Start switches moving the same amount - this ensures the same number
    %of samples are read each time, even though they vary by 1 or 2
    s_trim=[FreqStartSwitches+s_rem, FreqStartSwitches+s_int-s_rem];
    
    
else % we are ing singlefreqmode, and the interval is equal to the biggest gap between injection switches
    MultiFreqMode =0;
    s_int=max(max(diff(InjSwitches))); % interval in samples of each switch
    
    %calculate vector of trimmed segments - this is the gap between inj
    %switches
    s_trim=[InjSwitches(:,1)+s_rem, InjSwitches(:,1)+s_int-s_rem];
    
    
    % FIX LOADING DATA THAT ISNT THERE FOR VERY SHORT INJECTION
    s_trim (s_trim > length(Vdata)) = length(Vdata);
    
end




%% segment data into injections

%calculate segment width
seg_width=s_int-s_rem*2+1;

%create vector of voltages PRT x Voltage x CHN x Repeat
Vseg=nan(N_prt,seg_width,N_elec,N_rep);
Pseg=Vseg;

%% segment that shit

% iPrt=1;
rep=1;

if N_rep > 1
    
    for injcnt = 1:N_switch
        
%         if injcnt==N_switch
%            disp('blah'); 
%         end
        
        Vseg(iPrt,:,:,rep)=Vdata(s_trim(injcnt,1):s_trim(injcnt,2),:);
        Pseg(iPrt,:,:,rep)=Pdata(s_trim(injcnt,1):s_trim(injcnt,2),:);
        
        iPrt=iPrt+1; %update Prt pointer
        
        if iPrt > N_prt
            iPrt=1;
            rep=rep+1;
        end
    end
    
else
    %if there is only 1 repeat then ignore the 4th dimension
    for injcnt = 1:N_switch
        vtmp=Vdata(s_trim(injcnt,1):s_trim(injcnt,2),:);
        Vseg(iPrt,1:size(vtmp,1),:)=vtmp;
        
        ptmp=Pdata(s_trim(injcnt,1):s_trim(injcnt,2),:);
        Pseg(iPrt,1:size(ptmp,1),:)=ptmp;
        
        iPrt=iPrt+1; %update Prt pointer
        
        if iPrt > N_prt
            iPrt=1;
            rep=rep+1;
        end
        
        
    end
end

%report what the last protocol line was. it is incremented in the loop
%above so take one off it - KLUDGE
lastPrt=iPrt-1;

%if it is equal to zero that "wrap" it round to the last value
if lastPrt==0
    lastPrt=N_prt;
end

fprintf('Data segmented ok - last protocol line was %d of %d\n',lastPrt, N_prt);
end

