% PURPOSE:  converts ERPLAB data to a workable structure with the epp-tb
%           package.
%
% FORMAT
% ------
% study = epp_loaderplab(ALLERP)
% 
%
% INPUTS
% ------
% ALLERP    - structure built by ERPLAB package. All subject ERP sets
%             should be loaded before using this function. 
%
% OUTPUT
% ------
% study     - a structure orgenized by condition, with the following
%             feilds:
%                   1. Condition - name of condition (string)
%                   2. Data - ERP data, orgenized in 3D matrix
%                      (electrode*time*subject).
%                   3. timeLine - vector of time points corresponding to
%                      the time dimention in study.Data matrix.
%                   4. IDs - table containing IDs of subject in the
%                      condition.
%
% See also epp_loadegimat, epp_loadeeglab
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
06-03-2018  Removed samplingRate & baseLine fields from study struct.
25-11-2016  New function (written in MATLAB R2015a)
%}

function study = epp_loaderplab(ALLERP)

%% Import from ERPLABs ALLERP
ERPs.Condition = '';

for s = 1:length(ALLERP)
    fprintf('Loading file %d (of %d)...',s ,length(ALLERP))
    
    for c = 1:ALLERP(s).nbin
        if ~any(strcmpi(ALLERP(s).bindescr{c},{ERPs(:).Condition}))     % if condition hasn't been previously loaded:
            ERPs(end+1).Condition       = ALLERP(s).bindescr{c};
            ERPs(end).Data              = ALLERP(s).bindata(:,:,c);
            ERPs(end).timeLine          = ALLERP(s).times;
            ERPs(end).IDs               = table({ALLERP(s).filename},'VariableNames',{'ID'});
        else
            ind = strcmpi(ALLERP(s).bindescr{c},{ERPs(:).Condition});   % find where in the struct the condition was saved
            
            ERPs(ind).Data(:,:,end+1)	= ALLERP(s).bindata(:,:,c);     % append new data to end of data
            ERPs(ind).IDs{end+1,1}      = {ALLERP(s).filename};         % append new IDto end of ID list
        end
    end
    
    fprintf('Done!\n')
end

ERPs    = ERPs(2:end);  % remove first (blank) row.

study   = ERPs;         % save for export

fprintf('Done!\n')
end