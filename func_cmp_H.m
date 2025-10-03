%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the homogeneity statistic for the region 
% associated to a specific hindcast node and FS threshold 

% Hosking, J. R. M., & Wallis, J. R. (1997). Regional frequency analysis

% De Leo & Solari, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H1,H2,H3] = func_cmp_H(X)

%--------------------------------------------------------------------------
% FORMAT DATA FOR R PROCESSING
% 1) reshape matrix: each column indicates a site
Xt = X';

% 2) export temporary file
dlmwrite('temp_mt.dat',Xt,'Delimiter',' ','precision','%5.2f')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% RUN R ROUTINE
system('C:\Program Files\R\R-4.5.1\bin\x64\RScript run_03_HomTest.r')
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% IMPORT R RESULTS
D  = load('H_test.dat');
H1 = D(1);
H2 = D(2);
H3 = D(3);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% CLEAN WORKSPACE
system('del H_test.dat temp_mt.dat')
%--------------------------------------------------------------------------

return