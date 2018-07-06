% PURPOSE:  Retain subjects that have data in all spefified conditions.
%
%
% FORMAT
% ------
% [studyOut, nsubs] = epp_matchsubjects(studyIn,conditions)
%
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
%                 If left blank {}, subjects are matched across ALL
%                 conditions.
%
% See also epp_filter_by, epp_matchsubjects, epp_combineconds
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
06-07-2018  Rewrite to use epp_filter_by
09-06-2018  Added option to match across all conditions
13-05-2018  BUG FIX
05-03-2018  Added support for TF data
17-12-2017  Added output of number of subs
09-02-2017  Added support for more than two conditions
25-11-2016  New function (written in MATLAB R2015a)
%}

function [studyOut, nsubs] = epp_matchsubjects(studyIn,conditions)

if isempty(conditions)
    conditions = {studyIn.Condition};
end

%% Get only relevant conditions
cInd    = cellfun(@(x) find(strcmp(x,{studyIn(:).Condition})), conditions);
studyIn = studyIn(cInd);

%% Get only relevant subjects
innerjoined = studyIn(1).IDs(:,'ID');
for i = 2:length(conditions)
    innerjoined = innerjoin(innerjoined,studyIn(i).IDs(:,'ID'));
end

%% Filter

match_subs  = @(x) any(strcmp(x,innerjoined{:,'ID'}));
studyOut    = epp_filter_by(studyIn,'ID',match_subs,true); % use epp_filter_by
nsubs       = size(innerjoined,1);

end