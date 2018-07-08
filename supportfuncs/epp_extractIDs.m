% PURPOSE:  Extract data from study.IDs tables.
%
% FORMAT
% ------
% IDs_table = epp_extractIDs(study,conditions,varargin)
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
%                 If left blank, will return for all conditions.
%
% The available parameters are as follows:
%           'save'          - 'long' / 'wide'; will save the results in the
%                             current directory in the specified format.
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
07-07-2018  New function (written in MATLAB R2017a)
%}
function IDs_table = epp_extractIDs(study,conditions,varargin)

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addParameter(p,'save','no', @ischar);
parse(p ,study, conditions, varargin{:}); % validate

%% Get only relevant conditions

if isempty(conditions)
    conditions = {study.Condition};
end

cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study   = study(cInd);

%% combine tables

IDs_table = table();
for c = 1:length(study)
    study(c).IDs.Condition = repmat({study(c).Condition},1,size(study(c).IDs,1))';
    IDs_table = [IDs_table; study(c).IDs];
end
vNames = IDs_table.Properties.VariableNames;
vNames = setdiff(vNames,{'ID','Condition'});
IDs_table = IDs_table(:,{'ID','Condition',vNames{:}});

%% Save?
if any(strcmpi(p.Results.save, {'wide','long'}))
    fprintf('. Done!\n\n\n\nWriting to file..')
    % for each variable
    if strcmpi(p.Results.save,'wide')
        save_table = [];
        for v = 1:length(vNames)
            temp_t = IDs_table(:,{'ID','Condition',vNames{v}});
            temp_t = unstack(temp_t,vNames{v},'Condition');
            temp_t.Properties.VariableNames(2:end) = cellfun(@(x) [x '_' vNames{v}],temp_t.Properties.VariableNames(2:end),'UniformOutput',false);

            save_table = [save_table temp_t];
        end
    else
        save_table = IDs_table;
    end
    
    fn = ['ID_Data_' datestr(datetime, 'yyyymmdd_HHMMSS') '.csv']; % file name
    writetable(save_table,fn,'Delimiter',',','QuoteStrings',false) % write values
    
    fprintf('. Done!\n\n')
end
end