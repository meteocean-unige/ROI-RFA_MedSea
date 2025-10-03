%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs EVA following both RFA and at-site analysis
% for a target node in the Mediterranean Sea
% see:

% Regional Frequency Analysis of extreme waves based on
% Regions of Influence in the Mediterranean Sea

% De Leo & Solari, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

%% USER INPUT SETTINGS
% target node (RS: WGS84)
lonT = 18.569079;
latT = 35.766329;
 
nsim       = 1000;   % number of bootstrap for CI computation
alpha      = .05;    % CI level = (1-alpha)x100
pij_thresh = 0.65;   % cut-off for the detection of the ROI

Tr = [1.1:.1:2 3:10 10:10:100 200:100:1000]; % target ret. periods

%% LOAD DATA
load('clustered_peaks.mat')

lon   = swh_peaks.lon;    % longitude
lat   = swh_peaks.lat;    % latitude
Hs_AM = swh_peaks.AM;     % Hs AM series
ID    = swh_peaks.ID;     % hindcast ID
PK    = swh_peaks.peaks;  % peaks clustered through MeanShift++

%% RETURN LEVELS COMPUTATION

% select closest node to the target site
[~,irow] = min(pdist2([lonT latT],[lon lat]));

% ---------------------AT-SITE EVA ----------------------------------------
Ham_node = Hs_AM(irow,:);     % AM series at the target node
ny       = length(Ham_node);  % number of years

% empirical return periods (use Weibull plotting position)
Tr_hnd = 1./(1-[1:ny]./(ny+1));

% bootstap for CI computation
x_ci = nan(nsim,length(Tr));
for i = 1:nsim
    x                   = randsample(Ham_node,length(Ham_node),'true');
    [shMC,scMC,locMC]   = func_LMOM_GEV(x);
    x_ci(i,:)           = gevinv(1-1./Tr,shMC,scMC,locMC);
end
HsTr_CI_site = quantile(x_ci,[alpha/2 1-alpha/2],1);
% -------------------------------------------------------------------------

% ---------------------ROI/RFA EVA ----------------------------------------
% compute pij and detect ROI
PK(isnan(PK)) = 0;
PK(PK~=0)     = 1;
sumev  = PK + PK(irow,:);
isev   = sumev>0;
both   = sumev==2;
pij    = sum(both,2)./sum(isev,2);

id  = find(pij>pij_thresh);
id0 = find(id==irow);   % index of the target location within ROI's

Hs_AM_ROI = Hs_AM(id,:)';       % AM series within the ROI
pij_ROI   = pij(id,:);          % weights for regional moments computation
nroi      = size(Hs_AM_ROI,2);  % number of sites considered

% bootstap for CI computation accounting for inter-site correlation
idrnd = randi(ny,ny,nsim);
x_ci  = nan(nsim,length(Tr));
for idsim = 1:nsim
    amrnd             = Hs_AM_ROI(idrnd(:,idsim),:);
    [shMC,scMC,locMC] = func_LMOM_GEV_ROI(amrnd,pij_ROI,mean(amrnd(:,id0)));
    x_ci(idsim,:)     = gevinv(1-1./Tr,shMC,scMC,locMC);

end
HsTr_CI_roi = quantile(x_ci,[alpha/2 1-alpha/2],1);

%% PLOTS

% return level curves
fig=figure;
set(fig,'Position',[360 177 560 320]);

p1=semilogx(Tr,HsTr_CI_site,'-k');hold on
p2=semilogx(Tr,HsTr_CI_roi,'--r');hold on

scatter(Tr_hnd,sort(Ham_node),50,[1 1 1],'ko');

xlabel('T_R [y]')
ylabel('H_s [m]')

legend([p1(1) p2(1)],'at-site','ROI/RFA','Location','SouthEast')
legend boxoff

set(gca,'FontSize',12)

% ROI
fig=figure;
set(fig,'Position',[360   178   760   420])
geoscatter(lat(id),lon(id),5,pij(id),'Filled')
h = colorbar;
ylabel(h,'p_{ij}')
clim([pij_thresh 1]) % use caxis in older Matlab versions
geolimits([30 45], [-6 37])  
geobasemap grayland
set(gca,'FontSize',12)


