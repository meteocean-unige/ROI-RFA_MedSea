%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes regional L-moments of a dataset - these are next used 
% to compute GEV parameters (see Hosking & Wallis, 1997)

% De Leo & Solari, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kr,sr,mr,t,t3] = func_LMOM_GEV_ROI(data,w,l1o)

data = sort(data);
N    = size(data,1);
B0   = sum(data)./N;
B1   = sum(data(2:end,:).*(1:N-1)')./(N*(N-1));
B2   = sum(data(3:end,:).*(2:N-1)'.*(1:N-2)')./(N*(N-1)*(N-2));
B3   = sum(data(4:end,:).*(3:N-1)'.*(2:N-2)'.*(1:N-3)')./(N*(N-1)*(N-2)*(N-3));

l1 = B0;                        l1 = l1(:);
l2 = 2.*B1-B0;                  l2 = l2(:); t  = l2./l1;
l3 = 6.*B2-6.*B1+B0;            l3 = l3(:); t3 = l3./l2;
l4 = 20.*B3-30.*B2+12.*B1-B0;   l4 = l4(:); t4 = l4./l2;

% compute regional t statistics based on weighted means
tr  = sum(w.*t)./sum(w);
t3r = sum(w.*t3)./sum(w);

l2  = tr.*l1o;
C   = 2./(3+t3r)-log(2)/log(3);
kr  = 7.859.*C + 2.9554.*C.^2;

FUN = @(x) 2*(1-3^-x)/(1-2^-x)-3-t3r;
kr  = fzero(FUN,kr);

sr = l2.*kr./(1-2.^-kr)./gamma(1+kr);   % scale
mr = l1o - sr.*(1-gamma(1+kr))./kr;     % location
kr = -kr;                               % shape

return