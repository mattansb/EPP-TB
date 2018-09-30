% This function preforms a jackknife opperation on ERP data and is called
% internaly by measument supplumetary functions.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
14-04-2018  Function re-write for better preformance
19-12-2016  Add Jackknife adjustment as per: Smulders, F. T. (2010) 
25-11-2016  New function (written in MATLAB R2015a)
%}

function jmean = suppJackknife(mode,data,dim)
warning('suppJackknife is no longer supported. Use f_jackknife instead.')
jmean = f_jackknife(mode,data,dim);
end