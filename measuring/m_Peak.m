% PURPOSE:  measure latency \ amplitude (Peak)
%
% FORMAT
% ------
% [amp,lat] = m_Peak(data,timeWindow_ind,direction,local,times)
%
% See also epp_getLat, epp_getAmp
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
11-04-2018  Added help
07-04-2018  New function (written in MATLAB R2017a)
%}
function [amp,lat] = m_Peak(data,timeWindow_ind,direction,local,times)

try
    if direction < 0
        data = -1*data;
    end

    if local~=0 % if selected to find local peak
        ShiftRight  = data([end 1:(end-1)]);    % shift data to the right
        ShiftLeft   = data([2:end 1]);          % shift data to the left

        cLocal = data > ShiftRight & data > ShiftLeft; % find points which are GREATER than both adjecent points

        if local > 1 % so the same for the mean
            meanRight = tsmovavg(data, 's', local, 1);
            meanLeft = flip(tsmovavg(flip(data,1), 's', local, 1), 1);

            cLocal = data > meanRight & data > meanLeft & cLocal;
        end

        data(~cLocal) = nan; % remove point that do not mean crit

        timeWindow_ind = timeWindow_ind([local:(end-local)]); % cut down the time window
    end

    [amp, lat] = nanmax(data(timeWindow_ind)); % find amp and latency

    cut_times = times(timeWindow_ind);
    lat = cut_times(lat); % convert to time

    if direction < 0
        amp = -1*amp;
    end
    
    if isnan(amp)
        lat = nan;
    end
    
catch ME
    warning(ME.message)
    amp = nan;
    lat = nan;
end

end