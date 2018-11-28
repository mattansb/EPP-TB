% PURPOSE:  measure latency (Absolute criterion)
%
% FORMAT
% ------
% res = m_latAbsolute_criterion(data,timeWindow_ind,direction,times,criterion)
%
% See also epp_getLat
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
28-11-2018  Support for finding offset
11-04-2018  Added help
07-04-2018  New function (written in MATLAB R2017a)
%}
function res = m_latAbsolute_criterion(data,timeWindow_ind,direction,times,criterion,first_last)

try
    if direction < 0 % if max amp
        data        = -1*data;
        criterion   = abs(criterion);
    end
    
    cut_data    = data(timeWindow_ind);
    cut_times   = times(timeWindow_ind);
    
    ind = find(cut_data > criterion,1,first_last);   % find first time amp crosses criterion
    
    if ~isempty(ind)
        res = cut_times(ind);
    else
        res = nan;
    end
    
catch ME
    warning(ME.message)
    res = nan;
end

end