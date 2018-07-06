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

fn = fieldnames(studyIn);
has_erp     = any(strcmpi('Data',fn));
has_ersp    = any(strcmpi('ersp',fn));
has_itc     = any(strcmpi('itc',fn));

%% Get only relevant conditions
% cInd	= cellfun(@(x) find(ismember({studyIn(:).Condition}, x)), conditions);
cInd    = cellfun(@(x) find(strcmp(x,{studyIn(:).Condition})), conditions);
studyIn = studyIn(cInd);


%% Get only relevant subjects
for i = 2:length(conditions)
%     [~,ia,ib] = innerjoin(studyIn(1).IDs,studyIn(i).IDs);
    [~,ia,ib] = innerjoin(studyIn(1).IDs(:,1),studyIn(i).IDs(:,1));
    
    for j = 1:(i-1)
        studyIn(j).IDs  = studyIn(j).IDs(ia,:);
        
        if has_erp,  studyIn(j).Data = studyIn(j).Data(:,:,ia);   end
        if has_ersp, studyIn(j).ersp = studyIn(j).ersp(:,:,:,ia); end
        if has_itc,  studyIn(j).itc  = studyIn(j).itc(:,:,:,ia);  end
    end
    
    studyIn(i).IDs  = studyIn(i).IDs(ib,:);
    
    if has_erp,  studyIn(i).Data = studyIn(i).Data(:,:,ib);   end
    if has_ersp, studyIn(i).ersp = studyIn(i).ersp(:,:,:,ib); end
    if has_itc,  studyIn(i).itc  = studyIn(i).itc(:,:,:,ib);  end
end

studyOut = studyIn;

nsubs = size(studyOut(1).IDs,1);

end

