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


V=detrend(data);
N = length(V);


NFFT = max([2^16 2^nextpow2(length(V))]); % Next power of 2 from length of y
Y = fft(V,NFFT)/N;
F = Fs/2*linspace(0,1,NFFT/2+1);
Ymag=2*abs(Y(1:NFFT/2+1));


[~,maxw] = max(Ymag);
%find carrier one
Fc = F(maxw);
%display message to user






fprintf('Carrier frequency detected: Fc = %.2f Hz\r',Fc);

end

