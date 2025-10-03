%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function extracts SWH peaks based on the moving window approach
% De Leo & Solari, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dd    --> time array
% hh    --> initial series of SWH
% delta --> length of the time window sliding through the Hs series 

function [pk,thk,dk] = func_find_peaks(dd,hh,th,delta)

n = length(dd);

delta = delta*24;   % switch to hours to match hindicast data

% initialize arrays to speed up the computation
pk = nan(n,1); dk = nan(n,1); thk = nan(n,1);

k = 1;  % peaks counter

% loop through all data
for i = 1:(n-delta)
    
    endi  = i + delta;
    hdum  = hh(i:endi);
    hmax  = max(hdum);   

    imax  = i + find(hdum==hmax) - 1;    
    it    = round(i + (delta/2));
    
    % check if the peaks falls in the middle of the window - if so, store
    % the data in the respective arrays
    if it==imax(1) 

        pk(k)  = hmax;
        thk(k) = th(it);
        dk(k)  = datenum(dd(it));
        k = k+1;     

    end

end

pk(isnan(pk))   = [];
thk(isnan(thk)) = [];
dk(isnan(dk))   = [];

end


