% This function prepares ERP data to used in diffwave and LRP functions -
% and is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
09-06-2018  Added option to match across all conditions
13-05-2018  BUG FIX
05-03-2018  Added support for TF data
17-12-2017  Added output of number of subs
09-02-2017  Added support for more than two conditions
25-11-2016  New function (written in MATLAB R2015a)
%}

function [studyOut, nsubs] = suppMatchSubjects(studyIn,conditions)
warning('suppMatchSubjects is no longer supported. Use epp_matchsubjects instead.')
[studyOut, nsubs] = epp_matchsubjects(studyIn,conditions);

end

