% This function prepares ERP data to used in diffwave and LRP functions -
% and is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
% See also epp_getamplitude& & epp_getamplitude

%{
Change log:
-----------
17-12-2017  Added output of number of subs
09-02-2017  Added support for more than two conditions
25-11-2016  New function (written in MATLAB R2015a)
%}

function [studyOut, nsubs] = suppMatchSubjects(studyIn,conditions)

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
        studyIn(j).Data = studyIn(j).Data(:,:,ia);
    end
    
    studyIn(i).IDs  = studyIn(i).IDs(ib,:);
    studyIn(i).Data = studyIn(i).Data(:,:,ib);
end

studyOut = studyIn;

nsubs = size(studyOut(1).Data,3);

end

