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
% conditions    - cell list of conditions to be combined. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}),
%                 OR a regexp to match any of study(:).Condition.
%
% The available parameters are as follows:
%           'weighted'  - [defult: true] the combined waves are weighted
%                         according to nTrials variable in ID table.
%           'name'      - name of new combined condition. Defults is to
%                         concatenate the condition names.
%           'warning'   - {'on'|'off'} show warning when ID is missing from
%                         condition.
%
% See also epp_filter_by, epp_matchsubjects, epp_combineconds
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
06-07-2018  Add warning option.
29-05-2018  Add support for regexp + bug fix
20-03-2018  Bug fix.
05-03-2018  Added support for TF data
            Added ability to name output condition.
23-12-2017  New function (written in MATLAB R2017a)
%}
function out = epp_combineconds(study,conditions,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@(x) (iscellstr(x) & length(x)>=2) | ischar(x));
    addParameter(p,'weighted', true, @islogical)
    addParameter(p,'name', '', @ischar)
    addParameter(p,'warning', 'on', @ischar)
parse(p,study, conditions,varargin{:}); % validate

fn = fieldnames(study);
has_erp     = any(strcmpi('Data',fn));
has_ersp    = any(strcmpi('ersp',fn));
has_itc     = any(strcmpi('itc',fn));


%% Get only relevant conditions
if ischar(conditions)
    cInd        = cellfun(@(X) any(regexp(X,conditions)), {study.Condition});
    conditions  = {study(cInd).Condition};
    clear cInd
end

cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study   = study(cInd);

%% Subject List

IDs = {};
for c = 1:length(study)
    if ~any(isstr(study(c).IDs.ID))
        study(c).IDs.ID = cellfun(@(X) num2str(X),study(c).IDs.ID,'UniformOutput',false);
    end
    IDs = [IDs; study(c).IDs.ID];
end
try
    IDs = unique(IDs);
catch
    IDs = arrayfun(@(X) {X},unique(cell2mat(IDs)));
end

%% Combine
fprintf('Combining...\n')
for id = 1:length(IDs)
    % Make empty arrays
    if has_erp
        [nchans, ntime, ~] = size(study(1).Data);
        id_data = zeros([nchans, ntime, 1]);
    end

    if has_ersp
        [nchans, nfreq, ntime, ~] = size(study(1).ersp);
        id_ersp = zeros([nchans, nfreq, ntime, 1]);
    end

    if has_ersp
        [nchans, nfreq, ntime, ~] = size(study(1).itc);
        id_itc = zeros([nchans, nfreq, ntime, 1]);
    end
    
    % Combine conditions to arrays
    id_nTrials = [0];
    for c = 1:length(study)
        id_ind = find(strcmpi(IDs{id},study(c).IDs.ID));
        
        if isempty(id_ind)
            if strcmpi(p.Results.warning,'on')
                warning('ID %s missing data in condition %s',IDs{id},study(c).Condition);
            end
            continue
        end
        
        if p.Results.weighted
            w = study(c).IDs.nTrials(id_ind);
        else
            w = 1;
        end
        
        % append
        if has_erp,  id_data(:,:,end+1)    = study(c).Data(:,:,id_ind)*w;   end
        if has_ersp, id_ersp(:,:,:,end+1)  = study(c).ersp(:,:,:,id_ind)*w; end
        if has_itc,  id_itc(:,:,:,end+1)   = study(c).itc(:,:,:,id_ind)*w;  end
        
        id_nTrials(end+1)   = w;
        
        clear id_ind w
    end
    
    % Average for ID    
    nTrials(id) = sum(id_nTrials);
    if has_erp,  data(:,:,id)   = mean(id_data(:,:,2:end)/nTrials(id),3);   end
    if has_ersp, ersp(:,:,:,id) = mean(id_ersp(:,:,:,2:end)/nTrials(id),4); end
    if has_itc,  itc(:,:,:,id)  = mean(id_itc(:,:,:,2:end)/nTrials(id),4);  end
    
    clear id_data id_nTrials
end
fprintf('\b Done.\n')

%% Save

if any(cellfun(@(x) size(x,2),{study(:).IDs})>2)
    warning('Some ID data saved in IDs table has been lost')
end

out = study(1);
if isempty(p.Results.name)
    out.Condition = [conditions{:}];
else
    out.Condition = p.Results.name;
end

if has_erp,  out.Data = data; end
if has_ersp, out.ersp = ersp; end
if has_itc,  out.itc  = itc; end
out.IDs         = table(IDs,nTrials','VariableNames',{'ID' 'nTrials'});

end