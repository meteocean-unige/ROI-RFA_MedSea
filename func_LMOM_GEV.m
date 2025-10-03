%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code computes L-moments of a dataset - these are next used 
% to compute GEV parameters (see Hosking & Wallis, 1997)

% De Leo & Solari, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,SIG,MU] = func_LMOM_GEV(DATA)

DATA = reshape(DATA,length(DATA),1);
DATA = sort(DATA);
N    = size(DATA,1);

B0   = sum(DATA)./N;
B1   = sum(DATA(2:end,:).*(1:N-1)')./(N*(N-1));
B2   = sum(DATA(3:end,:).*(2:N-1)'.*(1:N-2)')./(N*(N-1)*(N-2));
B3   = sum(DATA(4:end,:).*(3:N-1)'.*(2:N-2)'.*(1:N-3)')./(N*(N-1)*(N-2)*(N-3));

LAM1 = B0;
LAM2 = 2.*B1-B0;
LAM3 = 6.*B2-6.*B1+B0;
LAM4 = 20.*B3-30.*B2+12.*B1-B0;

TAU  = LAM2./LAM1;
TAU3 = LAM3./LAM2;
TAU4 = LAM4./LAM2;

C    = 2./(3+TAU3)-log(2)/log(3);
K    = 7.859.*C + 2.9554.*C.^2;

for id = 1:numel(K)
    FUN   = @(x) 2*(1-3^-x)/(1-2^-x)-3-TAU3(id);
    K(id) = fzero(FUN,K(id));
end

SIG  = LAM2.*K./(1-2.^-K)./gamma(1+K); % scale
MU   = LAM1 - SIG.*(1-gamma(1+K))./K;  % location
K    = -K;                             % shape

return