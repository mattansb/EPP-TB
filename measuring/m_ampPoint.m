% PURPOSE:  measure amplitudes (point)
%
% FORMAT
% ------
% res = m_ampPoint(data,timeWindow_ind)
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
function res = m_ampPoint(data,timeWindow_ind)

try
    res = data(timeWindow_ind);
catch ME
    warning(ME.message)
    res = nan;
end

end