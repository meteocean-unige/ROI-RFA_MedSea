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

indir  = './archives/';
outdir = './archives/';

load([indir 'peaks_series_MS.mat']);
% It is assumed that the file ‘peaks_series_MS.mat’ contains the time 
% series of peaks for each node in the study area. The function 
% ‘func_find_peaks.m’ is included in the repository, which identifies the 
% series of peaks given a time series.

lon   = peaks.lon;
lat   = peaks.lat;
ID    = peaks.ID;
peaks = peaks.peaks; % This variable includes date and value of the peaks

% Annual Maxima at each node ----------------------------------------------
nnod = numel(peaks);
yys  = datevec(peaks{1}(:,1));
yys  = unique(yys(:,1));
yy0  = min(yys);
nyy  = numel(yys);
AM   = zeros(nnod,nyy);
tic
for idn = 1:nnod
    yynod = datevec(peaks{idn}(:,1));
    yynod = yynod(:,1);
    for idy = 1:nyy
        yyaux = yy0+idy-1;
        AM(idn,idy) = max(peaks{idn}(yynod==yyaux,2));
    end
end

%--------------------------------------------------------------------------
% Only the largest peaks are selected at each node for
% clustering, corresponding to an average of 12 peaks per year over the nyy 
% years of data.
nnod     = numel(peaks);
nall     = nnod*12*nyy;
peaksall = zeros(nall,3);
iddata2  = 0;
for idnod = 1:nnod
    iddata1 = iddata2+1;
    iddata2 = iddata2+504;
    [~,id] = sort(peaks{idnod}(:,2),'descend');
    peaksall(iddata1:iddata2,1:2) = peaks{idnod}(id(1:504),1:2);
    peaksall(iddata1:iddata2,3)   = idnod;
end

% Peaks are ordered by date -----------------------------------------------
[~,idsort] = sort(peaksall(:,1));
peaksall   = peaksall(idsort,:);

% Clustering algorithm: MeanShift++ ---------------------------------------
h   = [2 1 .4]; % Size of the cell in ºlon, ºlat and days.
b   = ones(5,5,5);
np  = size(peaksall,1);
y   = [lon(peaksall(:,3))-min(lon) ...
    lat(peaksall(:,3))-min(lat) ...
    peaksall(:,1)-min(peaksall(:,1))];
tol  = 1e-4;
maxd = 1;
iter = 0;
y0   = y;
while maxd>tol
    id  = round(y0./h)+1;
    s  = max(id);
    C  = zeros(s(1),s(2),s(3));
    S  = zeros([3,s]);
    for idp = 1:np
        C(id(idp,1),id(idp,2),id(idp,3)) = C(id(idp,1),id(idp,2),id(idp,3))+1;
        S(:,id(idp,1),id(idp,2),id(idp,3)) = ...
            S(:,id(idp,1),id(idp,2),id(idp,3)) + y0(idp,:)';
    end
    Cs = convn(C,b,'same');
    y1 = zeros(size(y0));
    for iddim = 1:3
        Ss          = convn(squeeze(S(iddim,:,:,:)),b,'same');
        ya          = Ss./Cs;
        y1(:,iddim) = ya(sub2ind(s,id(:,1),id(:,2),id(:,3)));
    end
    % -------
    % Checking for NaNs that might cause problems.
    idnan = find(isnan(y1));
    y1(idnan) = y0(idnan);
    % -------
    maxd = max(max(abs(y0-y1)));
    y0   = y1;
    iter = iter+1;
    disp([iter maxd])
end
[~,~,idc] = unique(y0,'rows');

% Builds a matrix of events (clusters) afecting each node -----------------
nevent = max(idc);
peaks = NaN(nnod,nevent);
for idevent = 1:nevent
    id = find(idc==idevent);
    peaks(peaksall(id,3),idevent) = peaksall(id,2);
end

% Prepares output ---------------------------------------------------------
swh_peaks.lon   = lon;
swh_peaks.lat   = lat;
swh_peaks.ID    = ID;
swh_peaks.AM    = AM;
swh_peaks.peaks = peaks;

% Saves output ------------------------------------------------------------
save([outdir 'clustered_peaks.mat'],'swh_peaks','-v7.3')