function [ compout ] = comparetank( simv,BVstruc,namestring,plotflag )
%comparetank makes the comparison plots for a collected dataset against
%simulated voltages
%   Currently assumes scousetom input

%% check inputs

if exist('namestring','var') == 0
    namestring='';
end

if exist('plotflag','var') ==0
    plotflag =[1 1 1];
end

%% Get mean experimental voltages

%this may have been calculated before
if ~(isfield(BVstruc,'BVave'))
    %if its the new cell type strucutre then take only the first freq *FOR
    %NOW*
    if iscell(BVstruc.BV)
      BVstruc.BVave=mean(BVstruc.BV{1},2);
    else
      BVstruc.BVave=mean(BVstruc.BV,2);
    end
end

%%

%convert collected data to mV
Vexp_full=(abs(BVstruc.BVave)/1e6).*1000;

%take abs of the simulated voltages in mV
Vsim_full=abs(simv).*1000;

%get keep_idx and rem
keep_idx=BVstruc.keep_idx;
rem_idx=BVstruc.rem_idx;

%get the electrode injections
injs=BVstruc.ExpSetup.Protocol;
chn=BVstruc.ExpSetup.Elec_num;
Ninj=size(injs,1);

%% Calculate the error between each measurement
e=abs(Vsim_full-Vexp_full);%./abs(Vsim_full);

%remove the error on injection channels 
e(rem_idx)=nan; %i had set this to zero because i didnt like how nans were handled by imagesc but that seems to be better now

%put this in a matrix of channel by electrode - from this we can see if
%there are patterns in the errors. Is one electrode consistantly measuring
%wrong? Or 
em=reshape(e,chn,Ninj);

%find injections for each elec
for ii=1:chn
    [tmp, ~]=find(injs == ii);
    elec_inj{ii}=tmp;
end

%find error for each injection
elec_e=zeros(Ninj,Ninj);

%find the mean of all the errors for each injection that includes each
%electrode - this will show if there are consistently bad *injections* 
for ii=1:chn
    tmp=(nanmean(em(:,elec_inj{ii})));
    elec_e(ii,1:length(tmp))=tmp;
end

%remove excess nans
elec_e=elec_e(:,any(elec_e));


%% percentage error
e_p=e(keep_idx)./Vsim_full(keep_idx)*100;


%% bin the errors within each mV as is easier to understand

%when checking the tank I found asking things like, does the error get
%bigger with smaller voltages easier to talk about when there was less than
%1000 points!

%take a bin of every mV
bins=[0:1:max(Vsim_full(keep_idx))];

%mean and std the data within each bin
for iBin=1:length(bins)-1
    e_pbin(iBin)=nanmean(e_p((Vsim_full(keep_idx) > (bins(iBin))) & (Vsim_full(keep_idx) < (bins(iBin+1)))));
    e_pbinstd(iBin)=nanstd(e_p((Vsim_full(keep_idx) > (bins(iBin))) & (Vsim_full(keep_idx) < (bins(iBin+1)))));
    e_bin(iBin)=nanmean(e((Vsim_full(keep_idx) > (bins(iBin))) & (Vsim_full(keep_idx) < (bins(iBin+1)))));
    e_binstd(iBin)=nanstd(e((Vsim_full(keep_idx) > (bins(iBin))) & (Vsim_full(keep_idx) < (bins(iBin+1)))));
end
%last bin is full of nonsense or is empty so ditch it
bins(end)=[];

%% Compute correlation coefficients for match

%the experimental and simulated should be the same. R matrix is 1 on diag
%as this is the autocorrelation, but the off diagonals should be as close
%to 1 if the data is good! Use the standard alpha of 95% for sig check, so
%if data is good, then the p matrix should be zeros on diags
[R,p]=corrcoef(Vsim_full(keep_idx),Vexp_full(keep_idx));



%% plot dat shit!

%first plot is the classic scatter plot
if plotflag(1) ==1
    figure;
    hold all
    plot(Vsim_full(keep_idx),Vexp_full(keep_idx),'o');daspect([1 1 1]);
    
    %this makes the ticks look nice by having a round number of them
    tickdiff=mode(diff(get(gca,'XTick')));
    extents=max([max(Vsim_full(keep_idx)) Vexp_full(keep_idx)']);
    extents=ceil((extents/tickdiff))*tickdiff;
    %plot the ideal y=x line
    plot([0 extents],[0 extents],'k:');
    hold off
    xlabel('Simulated mV')
    ylabel('Experimental mV')
    title([namestring ' - Comparison to simulation'])
    
    axis([0 extents 0 extents])
    
    % tightfig;
    
end

%the second plot is a bar of the mean of the errors for the injections on
%each electrode - this helps identify problem electrodes, which could then
%be related to the geometry of the tank
if plotflag(2) ==1
    
    for ii=1:size(elec_e,2);
        legstr{ii}=['Inj ' num2str(ii)];
    end
    
    figure;
    h=bar(elec_e,'stacked','EdgeColor','k','LineStyle','-');
    set(h,'EdgeColor','black')
    set(gca,'XTick',1:Ninj)
    xlim([0 chn+1])
    xlabel('Electrode')
    ylabel('Mean Error mV')
    title([namestring ' - errors on injections for each electrode'])
    legend(legstr);
    % tightfig;
end

%finally, plot the matrix of the errors, this should highlight consistnly
%bad electrodes or injections
if plotflag(3) ==1
    figure;
    imagesc(em)
    % caxis([0 20e-3])
    % caxis([0 0.0001])
    shading flat;
%     colormap((brewermap([],'*GnBu')))
    c=colorbar;
    
    
    %make the colorbar something normal like every 0.5mV, so it ends at a
    %round number
    tickdiff=mode(diff(get(c,'YTick')));
    extents=ceil((max(caxis)/tickdiff))*tickdiff;
    caxis([0 extents]);
    set(c,'YTick',0:tickdiff:extents);
    
    
    ylabel(c,'Error mV')
    title('Error Across Injections');
    xlabel('Injection');
    ylabel('Channel');
    % tightfig;
    
end
%% Make the output strucutre

compout.elec_e=elec_e;
compout.e=e;
compout.em=em;
compout.Vsim=Vsim_full;
compout.Vexp=Vexp_full;
compout.keep_idx=keep_idx;
compout.e_p=e_p;
compout.bins=bins;
compout.e_pbin=e_pbin;
compout.e_pbinstd=e_pbinstd;
compout.e_bin=e_bin;
compout.e_binstd=e_binstd;
compout.R=R;
compout.p=p;



end

