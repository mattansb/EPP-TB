% PURPOSE:  measure latency (Fractional area)
%
% FORMAT
% ------
% res = m_latFractional_area(data,timeWindow_ind,direction,boundary,times,percentage)
%
% See also epp_getLat
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
02-01-2019  (re)Added linear interpolation
11-04-2018  Added help
07-04-2018  New function (written in MATLAB R2017a)
%}
function res = m_latFractional_area(data,timeWindow_ind,direction,boundary,times,percentage)

try
    data = data-boundary;

    % Rectify area
    if direction == 0 % all area
        data = abs(data); % make all amps positive!
    elseif direction > 0 % only positive area
        data(data < 0) = 0;
    elseif direction < 0 % only negative area
        data(data > 0)  = 0;
        data            = -1*data;
    end

    % turn amps in time window to cumulative area: 
    cut_data    = data(timeWindow_ind);
    cut_data    = cumsum(cut_data);
    cut_data    = cut_data/cut_data(end);
    cut_times   = times(timeWindow_ind);

    % compute criterion = %*area
    ind = find(cut_data > percentage,1);   % find first time amp crosses criterion
    ind = [ind-1,ind];
    
    if ~isempty(ind)
        res = interp1(cut_data(ind),cut_times(ind),percentage); % linear interpolation
%         res = cut_times(ind);
    else
        res = nan;
    end
    
catch ME
    warning(ME.message)
    res = nan;
end

end