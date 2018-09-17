function [ prt_full,keep_idx,rem_idx,Elec_inj ] = ScouseTom_data_findprt( InjectionPairs,N_elec )
%[prt_full,keep_idx,rem_idx,Elec_inj] = ScouseTom_data_findprt( InjectionPairs,N_elec )
%SCOUSETOM_DATA_FINDPRT Creates full measurement protocol, with rem_idx and
%keep_idx from injection pairs and number of electrodes. Also provides the
%channels which represent the Injection channels for estimating the
%impedance


%% Create full protocol
vp=(1:N_elec)';% positive voltage channel 1  to Number of electrodes
vm=ones(size(vp))*(N_elec+1); %always against a ground electrode

%make the entire protocol each line is INJ+ INJ- MEAS+ MEAS-
prt_mat=[];
for iii=1:size(InjectionPairs,1)
    temp=[repmat(InjectionPairs(iii,:),N_elec,1) vp vm];
    prt_mat=[prt_mat ; temp];
end
%% keep_idx and rem_idx
%find remove index, any protocol lines including the injetion channels are
%"bad"
prt_full=prt_mat;
prt=prt_full;
rem_idx=[];
for iPrt = 1:size(prt,1)
    if any(ismember(prt_full(iPrt,1:2),prt(iPrt,3:4))) ==1
        rem_idx=[rem_idx,iPrt];
    end
end
%keep index is anything that we *dont* remove
keep_idx=setdiff(1:length(prt_full),rem_idx);

N_prt=size(InjectionPairs,1);

%% Find instances of injecting on each electrode
%Electrode Injections
Elec_inj_mat=nan(N_elec,N_prt);

for iPrt = 1:N_prt
    Prt_cur=InjectionPairs(iPrt,:);
    start_idx=((iPrt-1)*N_elec);
    BV_chn=start_idx+Prt_cur;
    Elec_inj_mat(Prt_cur,iPrt)=BV_chn;
end

Elec_inj_mat=sort(Elec_inj_mat,2);

for iElec=1:N_elec
    Elec_inj{iElec}=Elec_inj_mat(iElec,~isnan(Elec_inj_mat(iElec,:)));
end
end