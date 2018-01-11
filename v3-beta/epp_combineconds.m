% PURPOSE:  Combine conditions into single condition.
%
% 
% FORMAT
% ------
% out = epp_combineconds(study,conditions,varargin)
%
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
%
% The available parameters are as follows:
%           'weighted'  - [defult: true] the combined waves are weighted
%                         according to nTrials.
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
2DO
----
make also for ERSP and ITC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Change log:
-----------
23-12-2017  New function (written in MATLAB R2017a)
%}
function out = epp_combineconds(study,conditions,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addParameter(p,'weighted', true, @islogical)
parse(p,study, conditions,varargin{:}); % validate

%% Get only relevant conditions
cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study   = study(cInd);

%% Subject List

IDs = {};
for c = 1:length(study)
    IDs = {IDs{:}, study(c).IDs.ID{:}};
end
IDs = unique(IDs);

%% Combine
[nchans, ntime, ~] = size(study(1).Data);

for id = 1:length(IDs)
    id_data = zeros([nchans, ntime, 1]);
    id_nTrials = [0];
    for c = 1:length(study)
        id_ind = find(strcmpi(IDs{id},study(c).IDs.ID));
        
        if isempty(id_ind)
            msg = fprintf('ID %s missing data in condition %s',IDs{id},study(c).Condition);
            warning(msg)
            continue
        end
        
        if p.Results.weighted
            w = study(c).IDs.nTrials{id_ind};
        else
            w = 1;
        end
        
        id_data(:,:,end+1) = study(1).Data(:,:,id_ind)*w;
        id_nTrials(end+1) = w;
        
        
        clear id_ind w
    end
    data(:,:,id) = mean(id_data(:,:,2:end),3);
    nTrials(id) = sum(id_nTrials);
    clear id_data id_nTrials
end

%% Save

if any(cellfun(@(x) size(x,2),{study(:).IDs})>2)
    warning('Some ID data saved in IDs table has been lost')
end

out             = study(1);
out.Condition   = [conditions{:}];
out.Data        = data;
out.IDs         = table(IDs',nTrials','VariableNames',{'ID' 'nTrials'});

end