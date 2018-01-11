% PURPOSE:  creates an ERP structured containing a DiffWave conditions. The
%           function first makes sure both conditions have the same
%           subjects (if doesnt, only data from matching subjects is used).
%
% 
% FORMAT
% ------
% DIFF = epp_diffwave(study, conditions)
%
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
22-03-2017  Better naming of diffwave condition
25-11-2016  New function (written in MATLAB R2015a)
%}

function DIFF = epp_diffwave(study, conditions)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@(x) iscellstr(x) && length(x)==2);
parse(p,study, conditions); % validate

%% Prepare Data

DIFF = suppMatchSubjects(study,conditions);

%% Calculate Diffwave

% DIFF(1).Condition   = strjoin(conditions,'-');
DIFF(1).Condition   = [conditions{1} '_' conditions{2} '_diff'];
DIFF(1).Data        = DIFF(1).Data-DIFF(2).Data;

DIFF = DIFF(1);

end