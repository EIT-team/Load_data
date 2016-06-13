%% This function examines EIT data and determines which channels are being used for injection
% and what frequency is being injected. This is acheived by finding the
% most prominent frequency at each electrode and then assuming that the two
% largest electrodes for a given frequency are the injection electrodes.
% Intended for use for parallel injection, but should also work for single,
% but not as useful in that case.

%Inputs - Data: values from EEG amplifier (BioSemi, Actichamp) etc, with
%N_elecs * N_Samples size.
%       Fs: Sampling frequency of EEG amp

% Tom Dowrick Nov 2015

function [Prt, Freqs] = Find_Injection_Freqs_And_Elecs(Data,Fs)

N_elecs = size(Data,2);
%Stops annoying error message spamming the screen
warning('off','MATLAB:colon:nonIntegerIndex')

%Get most prominent frequency on each electrode
for i = 1:N_elecs
    [Pxx,w] = pwelch(Data(:,i),Fs,.95,[],Fs);
    w_ind = find(w>105);
    [a,maxw] = max(Pxx(w_ind));
    
    
    Elec_Freq_Mag(i,:) = [i w(w_ind(maxw)) a]; %Stores electrode number, frequency, and magnitude at that frequency
    clear Pxx w w_ind max
end

%Get unique frequencies, should correspond to measurement frequencies used!
Freqs = unique(Elec_Freq_Mag(:,2));
N_Freqs = length(Freqs);

%Calculate the injection electodes for each frequency, by looking for the
%channels with the two largest magnitdues for a given frequency.

for i =1:N_Freqs
    
    A = Elec_Freq_Mag(Elec_Freq_Mag(:,2) == Freqs(i),:); %Get only the values for Freq(i)
    [~,ind_max] = max(A(:,3)) ;     %index with maximum amplitude
    Prt(i,1) = A(ind_max,1)  ;       %Store electrode max value corresponds to
    A(ind_max,3) = -Inf;            %Need to find 2nd largest value, set previous max to some low value
    [~,ind_max] = max(A(:,3)) ;      %Find 2nd biggest
    Prt(i,2) = A(ind_max,1) ;        %Store 2nd electrode
    
end

for i = 1:size(Freqs)
    disp(['Detected ' num2str(Freqs(i)) 'Hz between electrodes ' num2str(Prt(i,1)) ' and ' num2str(Prt(i,2))])
end

end