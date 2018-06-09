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
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be subtracted (second will be
%                 subtracted from the first). Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
%
% The available parameters are as follows:
%           'name'      - name of new combined condition. Defults is to
%                         concatenate the condition names.
%
% See also epp_combineconds, epp_LRP, epp_GFP, epp_makegrands
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
20-03-2018  Fix documentation.
05-03-2018  Added support for TF data
            Added ability to name output condition.
22-03-2017  Better naming of diffwave condition
25-11-2016  New function (written in MATLAB R2015a)
%}

function DIFF = epp_diffwave(study, conditions,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@(x) iscellstr(x) && length(x)==2);
    addParameter(p,'name', '', @ischar)
parse(p,study, conditions,varargin{:}); % validate

%% Prepare Data

DIFF = epp_matchsubjects(study,conditions);

fn = fieldnames(study);
has_erp     = any(strcmpi('Data',fn));
has_ersp    = any(strcmpi('ersp',fn));
has_itc     = any(strcmpi('itc',fn));

%% Calculate Diffwave

% Condition name
if isempty(p.Results.name)
    DIFF(1).Condition = [conditions{1} '_' conditions{2} '_diff'];
else
    DIFF(1).Condition = p.Results.name;
end

% Prep data:
if has_erp,  DIFF(1).Data = DIFF(1).Data - DIFF(2).Data; end
if has_ersp, DIFF(1).ersp = DIFF(1).ersp - DIFF(2).ersp; end
if has_itc,  DIFF(1).itc = DIFF(1).itc   - DIFF(2).itc; end

DIFF(1).IDs = DIFF(1).IDs(:,1);

DIFF = DIFF(1);

end