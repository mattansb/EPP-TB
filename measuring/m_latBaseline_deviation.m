% PURPOSE:  measure latency (Baseline deviation)
%
% FORMAT
% ------
% res = m_latBaseline_deviation(data,timeWindow_ind,direction,baselineOS,times,SDcriterion,first_last)
%
% See also epp_getLat
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
02-01-2019  (re)Added linear interpolation
28-11-2018  Support for finding offset
11-04-2018  Added help
07-04-2018  New function (written in MATLAB R2017a)
%}
function res = m_latBaseline_deviation(data,timeWindow_ind,direction,baselineOS,times,SDcriterion,first_last)

try
    if direction < 0 % if max amp
        data = -1*data;
    end
    
    if baselineOS < 0
        baseline_ind = times >= baselineOS & times < 0;
    elseif baselineOS > 0
        baseline_ind = times <= baselineOS & times > 0;
    end

    Cr          = std(data(baseline_ind))*SDcriterion; % compute criterion = n*SD     
    cut_data    = data(timeWindow_ind);
    cut_times   = times(timeWindow_ind);
    
    ind = find(cut_data > Cr,1,first_last);   % find first / last time amp crosses criterion
    
    if ~isempty(ind)
        switch first_last
            case 'first'
                ind = [ind-1,ind];
            case 'last'
                ind = [ind,ind+1];
        end
        res = interp1(cut_data(ind),cut_times(ind),Cr); % linear interpolation
%         res = cut_times(ind);
    else
        res = nan;
    end
    
catch ME
    warning(ME.message)
    res = nan;
end

end