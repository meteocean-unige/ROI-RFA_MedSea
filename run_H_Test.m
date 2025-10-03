%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code checks on the homogeneity of a ROI for a target node
% depending on the threshold values of the FS index
% reference is made to nodes used in the article (see Table 2)

% it uses func_cmp_H, that requires to have either R or RStudio installed 
% and the associated R script place in the same directory-otherwise add the
% path to the script

% See:

% Regional Frequency Analysis of extreme waves based on
% Regions of Influence in the Mediterranean Sea

% De Leo & Solari, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

% pij cut-off for the detection of the ROI 
pij_thresh = 0.1:0.1:0.9; 

% target node (RS: WGS84)
lonT = 18.569079;
latT = 35.766329;

% data for all Mediterranean Sea
load('clustered_peaks.mat')

lon   = swh_peaks.lon;    % longitude
lat   = swh_peaks.lat;    % latitude
Hs_AM = swh_peaks.AM;     % Hs AM series
ID    = swh_peaks.ID;     % hindcast ID
PK    = swh_peaks.peaks;  % clustered peaks

% build 0/1 series
PK(isnan(PK)) = 0;
PK(PK~=0)     = 1;

%----------------------------------------------------------------------
% detect reference node
[~,irow] = min(pdist2([lonT latT],[lon lat]));

lon_NODE = lon(irow);
lat_NODE = lat(irow);
ID_NODE  = ID(irow);

%----------------------------------------------------------------------
% compute Weiss index
sumev = PK + PK(irow,:);
isev  = sumev>0;
both  = sumev==2;
pij = sum(both,2)./sum(isev,2);

%--------------------------------------------------------------------
% LOOP THROUGH THE PIJ THRESHOLDS
H   = nan(length(pij_thresh),3);  % heterogeneity measure for the region
nn  = nan(length(pij_thresh),1);  % number of nodes
for ik = 1:length(FS_thresh)
    
    % computation of the H statistic
    % retain only nodes matching the FS condition
    Hs_AMC = Hs_AM(dij>pij_thresh(ik),:);
    [H1,H2,H3] = func_cmp_H(Hs_AMC);
    
    H(ik,1) = H1;
    H(ik,2) = H2;
    H(ik,3) = H3;
    
    nn(ik)  = size(Hs_AMC,1);
end

% Plot p-H
figure
plot(pij_thresh,H,'o-');
xlabel('p'); ylabel('H'); legend('H1','H2','H3');

% Plot p-nn
figure
plot(pij_thresh,nn,'o-');
xlabel('p'); ylabel('number of node in region');