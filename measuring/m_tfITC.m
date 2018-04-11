% PURPOSE:  measure TF (mean ITC)
%
% FORMAT
% ------
% res = m_tfITC(data,timeWindow_ind)
%
% See also epp_getTF
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
11-04-2018  Added help
07-04-2018  New function (written in MATLAB R2017a)
%}
function res = m_tfITC(data,timeWindow_ind)

try
    res = mean(abs(data(timeWindow_ind))); % this is the itc data we want    
catch ME
    warning(ME.message)
    res = nan;
end

end