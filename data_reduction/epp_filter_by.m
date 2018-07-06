% PURPOSE:  Filter data in study by some variable in study.IDs.
%
% 
% FORMAT
% ------
% out = epp_filter_by(study,varname,anon_fun,should_keep)
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% varname       - (string) name of variable in study.IDs table.
% anon_fun      - an anonymous function function. This function will be
%                 applied to the column specifyed by `varname`. It should
%                 return a logical value. 
% should_keep   - (logical). Does the anonymous function indicate rows to
%                 keep [true]? Or to remove [false]? 
%
%
% See also epp_filter_by, epp_matchsubjects, epp_combineconds
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
06-07-2018  New function (written in MATLAB R2017a)
%}
function study = epp_filter_by(study,varname,anon_fun,should_keep)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'varname',@ischar);
    addRequired(p,'anon_fun',@(x) isa(x,'function_handle'));
    addRequired(p,'should_keep',@islogical);
parse(p,study, varname,anon_fun,should_keep); % validate

fn = fieldnames(study);
has_erp     = any(strcmpi('Data',fn));
has_ersp    = any(strcmpi('ersp',fn));
has_itc     = any(strcmpi('itc',fn));

%% Filter
for c = 1:length(study)
    inds = arrayfun(anon_fun,study(c).IDs{:,varname});
    
    if ~should_keep, inds = ~inds; end
    
    % Remove from ID table
    study(c).IDs = study(c).IDs(inds,:);
    
    % Remove from data
    if has_erp,     study(c).Data   = study(c).Data(:,:,inds);      end
    if has_ersp,    study(c).ersp   = study(c).ersp(:,:,:,inds);    end
    if has_itc,     study(c).itc    = study(c).itc(:,:,:,inds);     end
end

end