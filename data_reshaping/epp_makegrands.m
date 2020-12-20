% PURPOSE:  Collapes all data across subjects. This is useful for plotting
%           data of large data sets.
%
%
% FORMAT
% ------
% study = epp_makegrands(study)
%
%
% See also epp_LRP, epp_GFP, epp_makegrands, epp_diffwave
%
%
% Author: Mattan S. Ben Shachar & Rachel Rac, BGU, Israel

%{
Change log:
-----------
23-04-2018  Minor fix
16-04-2018  New function (written in MATLAB R2017a)
%}
function study = epp_makegrands(study)
for c = 1:length(study)
    %% IDs
    study(c).IDs = table(000,size(study(c).IDs,1),'VariableNames',{'ID' 'nTrials'});
    
    %% Data
    if isfield(study,'Data')
        study(c).Data = mean(study(c).Data,3);
    end

    if isfield(study,'ersp')
        study(c).ersp = mean(study(c).ersp,4);
    end

    if isfield(study,'itc')
        study(c).itc = mean(study(c).itc,4);
    end
end