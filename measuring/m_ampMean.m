% PURPOSE:  measure amplitudes (mean)
%
% FORMAT
% ------
% res = m_ampMean(data,timeWindow_ind)
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
function res = m_ampMean(data,timeWindow_ind)

try
    res = mean(data(timeWindow_ind)); % compute mean
catch ME
    warning(ME.message)
    res = nan;
end

end