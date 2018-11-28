% PURPOSE:  measure latency (Relative criterion)
%
% FORMAT
% ------
% res = m_latRelative_criterion(data,timeWindow_ind,direction,local,times,percentage)
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
function res = m_latRelative_criterion(data,timeWindow_ind,direction,local,times,percentage,first_last)

try
    [amp,lat] = m_Peak(data,timeWindow_ind,direction,local,times);
    
    if isnan(amp)
        res = nan;
        return
    end
    
    if direction < 0 % if max amp
        data    = -1*data;
        amp     = -1*amp;
    end
    
    Cr              = amp*percentage;       % compute criterion = %*peak_amp
    lat             = find(times==lat,1);   % convert latency back to index
    
    switch first_last
        case 'first'
            data(lat:end) = nan; % remove later latencies
        case 'last'
            data(1:lat) = nan; % remove earlier latencies
    end
    
    ind = find(data > Cr,1,first_last); % find where amp is smaller \ larger than Criterion
    
    
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