function [ Fc ] = ScouseTom_data_GetCarrier( data,Fs )
%get carrier frequency of data
% Assumes only one carrier frequency 
%   finds the largest frequency bin in a data set. Zero pads the result to
%   give smaller frequency bins for short signals for a more accurate
%   result.
% Adapted from G Dragons Code by the sonorous and majestic Jimmy

% this is not a very robust bit of code at the moment, if data is fucked
% then this doesnt know, and at the moment there are no checks to see if it
% matches the *expected* freq




%find spectrum of signal with no overlaps or windowing
    [Pxx,w] = pwelch(detrend(data),[],0,2^16,Fs);
    % we only want frequencies above 3 Hz 
    w_ind = find(w>3);
    %find the biggest one
    [~,maxw] = max(Pxx(w_ind));
    %find carrier one
    Fc = w(w_ind(maxw));
    %display message to user
   fprintf('Carrier frequency detected: Fc = %.2f Hz\r',Fc);    
    
end

