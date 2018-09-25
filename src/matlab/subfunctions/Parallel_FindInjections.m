function [Prt, Freqs] = Parallel_FindInjections(Data,Fs,Chn_labels)
% [Prt, Freqs] = Parallel_FindInjections(Data,Fs,Chn_labels)
% This function examines EIT data and determines which channels are being used for injection
% and what frequency is being injected. This is acheived by finding the
% most prominent frequency at each electrode and then assuming that the two
% largest electrodes for a given frequency are the injection electrodes.
% Intended for use for parallel injection, but should also work for single,
% but not as useful in that case.
%
%Inputs - Data: values from EEG amplifier (BioSemi, Actichamp) etc, with
%N_elecs * N_Samples size.
%       Fs: Sampling frequency of EEG amp
%       Chn_labels: the correct electrode numbering, from HDR.Label. This
%       is so the channels are called the same as the *physical* channels
%       on the actichamp. Used for output text only 
%
%%
max_length = 2^20; %mamximum length of data to use, must be power of 2
Data=detrend(Data); % remove DC offset

if length(Data) > max_length
    Data=Data(1:max_length,:);
end

N = length(Data);

NFFT = max([max_length 2^nextpow2(length(Data))]); % Next power of 2 from length of y
N_elecs = size(Data,2);

if exist('Chn_labels','var') ==0 || isempty(Chn_labels)
    Chn_labels =1:N_elecs;
end
%%
Freq_Idx=nan(N_elecs,2);

%Get most prominent frequency on each electrode
for ichn = 1:N_elecs
    Y = fft(Data(:,ichn),NFFT)/N; % corrected power
    F = Fs/2*linspace(0,1,NFFT/2+1); %one side of frequency
    
    Ymag=2*abs(Y(1:NFFT/2+1)); % magnitude
    
    [~,maxw] = max(Ymag);
    %find carrier one
    Fc = F(maxw);
    Fc = round(Fc,2);
    
    Freq_Idx(ichn,:) = [Fc maxw]; %Stores electrode number, frequency, and magnitude at that frequency
end

%Get unique frequencies, should correspond to measurement frequencies used!
Freqs = unique(Freq_Idx(:,1));
Idx_F = unique(Freq_Idx(:,2));
N_freqs = length(Freqs);

%% Find the magnitude for each frequency on each elec

% now we know what the unique injection frequencies are, find the amplitude
% of all of these frequencies on all of the channel. Then take the biggest
% two channels per freq

Freq_mag=nan(N_elecs,N_freqs);
for ichn = 1:N_elecs
    Y = fft(Data(:,ichn),NFFT)/N; % corrected power
    Ymag=2*abs(Y(1:NFFT/2+1)); % magnitude
    
    for iFreq = 1:N_freqs
        Freq_mag(ichn,iFreq)=Ymag(Idx_F(iFreq));
    end
    
end

% find biggest two channels for each freq
for iFreq = 1:N_freqs
    
    curMag=Freq_mag(:,iFreq);
    
    [M,idx]=sort(curMag,1,'descend');
    
    Prt(iFreq,:)=idx(1:2);
    
end
% make sure the protocols are in order
Prt=sort(Prt,2);
%output to user
for ichn = 1:size(Freqs)
    disp(['Detected ' num2str(round(Freqs(ichn))) ' (' num2str(round(Freqs(ichn),-2)) ') Hz between electrodes ' num2str(Chn_labels(Prt(ichn,1))) ' and ' num2str(Chn_labels(Prt(ichn,2)))])
end

end