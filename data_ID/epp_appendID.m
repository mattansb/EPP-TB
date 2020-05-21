% PURPOSE:  Join ID data from a table
%
% FORMAT
% ------
% study = epp_appendID(study, T, varargin)
% 
%
% INPUTS
% ------
% study     - structure.
% T         - A table to be left joined to study.IDs.
%
% The available parameters are as follows:
%           'Keys'      - An index or name of the column by which to join.
%                         Default is 'ID'.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
%{
Change log:
-----------
08-05-2020  New function (written in MATLAB R2017b)
%}
function study = epp_appendID(study, T, varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'T',@istable);
    
    addParameter(p,'Keys','ID');
parse(p, study, T, varargin{:}); % validate

%% Left joint

for c = 1:length(study)
    in_table = study(c).IDs;
    
    out_table = outerjoin(study(c).IDs, T, 'type', 'left',...
        'MergeKeys', true,'Keys', p.Results.Keys);
    
    if height(out_table) ~= height(in_table)
        fprintf('\nIn table (%s):\n', study(c).Condition)
        disp(in_table)
        fprintf('\nOut table:\n')
        disp(out_table)
        error('Something went wrong - the out-table has a different number of rows than the in-table.')
    end
    
    study(c).IDs = out_table;    
end

end