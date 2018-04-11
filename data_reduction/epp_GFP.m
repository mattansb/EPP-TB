% PURPOSE:  creates an ERP structured containing data that has been reduced
%           across electrodes, in one of two methods:
%               + [defult] Global Field Power (GFP)
%                 as per Lehmann & Skrandies (1984).
%               + Scalp Effect Strength (SES)
%                 as per  Milivojevic, Johnson,Hamm, & Corballis (2003).
%
%
% 
% FORMAT
% ------
% GFP = epp_GFP(study,conditions,varargin)
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
%           'SES'     - if a string is given, the Scalp Effect Strength
%                       (SES) is calculated, as per  Milivojevic, Johnson,
%                       Hamm, & Corballis (2003).
%
% See also epp_combineconds, epp_LRP, epp_diffwave
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
21-12-2017  New function (written in MATLAB R2017a)
%}
function GFPs = epp_GFP(study,conditions,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addParameter(p,'SES', 'off', @ischar)
parse(p,study,conditions,varargin{:}); % validate

switch p.Results.SES
    %% Get GFP
    case 'off'
        cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
        study   = study(cInd);

        nChan = size(study(1).Data,1);
        for c = 1:length(study)
            temp_data = study(c).Data;
            for i = 1:nChan
                temp_GFP(:,:,:,i) = sum((temp_data(i,:,:)-temp_data(:,:,:)).^2,1);
            end
            temp_GFP = sqrt(sum(temp_GFP,4)*(1/(2*nChan)));
            if ndims(temp_GFP)==2
                temp_GFP = temp_GFP';
            end
            study(c).Data = temp_GFP;
            clear temp_GFP temp_data
        end
        GFPs = study;
    otherwise
        for c = 1:length(conditions)
            GFPs(c) = epp_diffwave(study,{p.Results.SES,conditions{c}});
            GFPs(c).Condition = conditions{c};
            temp_data = GFPs(c).Data;
            temp_data = sum(temp_data.^2,1);
            GFPs(c).Data = temp_data;
        end
end


end