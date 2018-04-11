% PURPOSE:  measure amplitudes (integral)
%
% FORMAT
% ------
% res = m_ampIntegral(data,timeWindow_ind,samplingRate)
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
function res = m_ampIntegral(data,timeWindow_ind,samplingRate)

try
    res = sum(data(timeWindow_ind))/samplingRate; % compute integral mv/sec
catch ME
    warning(ME.message)
    res = nan;
end

end