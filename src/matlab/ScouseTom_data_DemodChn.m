function [ Vdemod,Pdemod ] = ScouseTom_data_DemodChn( V,B,A )
%ScouseTom_data_DemodChn Demodulates each channel one at a time
%   Detailed explanation goes here

%% data info

 N_elec=size(V,2);
% N_prt=size(Vseg,1);
% N_rep=size(Vseg,4);


%% demodulate data

Vdemod=nan(size(V));
Pdemod=nan(size(Vdemod));

    for iElec=1:N_elec

            [Vdemod(:,iElec), Pdemod(:,iElec)]=ScouseTom_data_DemodHilbert(V(:,iElec),B,A);

    end
% disp('Demodulation done');


end

