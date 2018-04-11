% PURPOSE:  measure amplitudes (area)
%
% FORMAT
% ------
% res = m_ampArea(data,timeWindow_ind,direction,samplingRate,boundary)
%
% See also epp_getAmp
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
11-04-2018  Added help
07-04-2018  New function (written in MATLAB R2017a)
%}
function res = m_ampArea(data,timeWindow_ind,direction,samplingRate,boundary)

try
    data = data - boundary;
    if direction > 0      % only use positive values
        data(data < 0) = 0;
    elseif direction < 0  % only use negative values
        data(data > 0) = 0;
    end

    res = sum(abs(data(timeWindow_ind)))/samplingRate; % compute area mv/sec
catch ME
    warning(ME.message)
    res = nan;
end



end