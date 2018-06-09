% PURPOSE:  creates an ERP structured containing a DiffWave conditions. The
%           function first makes sure both conditions have the same
%           subjects (if doesnt, only data from matching subjects is used).
%
% 
% FORMAT
% ------
% LRP = epp_LRP(study, conditions, electrodes)
%
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - two options:
%                   1. only ONE condition - LRP is produced from only one
%                      condition. 
%                   2. TWO conditions - LRP is produced by avereging both
%                      conditions' (opposite) LRPs.
% electrodes    - 2*E matrix of electrodes, with each pair to be subtracted
%                 from one another (e.g. [87 85, 92; 65, 55 ,54]).
%
%
% OUTPUT
% ------
% NOTE that output has only the number of selected electrodes (ordered as
% they were input).
%
% See also epp_combineconds, epp_GFP, epp_diffwave, epp_makegrands
% 
% 
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
25-11-2016  New function (written in MATLAB R2015a)
%}

function LRP = epp_LRP(study, conditions, electrodes)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@(x) iscellstr(x) && length(x)>=2); % <<<<<<<<<<<<<<<<?
    addRequired(p,'electrodes',@isnumeric);
parse(p,study, conditions, electrodes); % validate

%% Make LRP

if length(conditions)==1
    LRP                 = study;
    LRP.Condition       = [LRP.Condition '_LRP'];
    LRP.Data            = LRP.Data([electrodes(1,:)],:,:)-LRP.Data([electrodes(2,:)],:,:);
else
    LRP                 = epp_matchsubjects(study,conditions);
    LRP(1).Data         = LRP(1).Data([electrodes(1,:)],:,:)-LRP(1).Data([electrodes(2,:)],:,:);
    LRP(2).Data         = LRP(2).Data([electrodes(2,:)],:,:)-LRP(2).Data([electrodes(1,:)],:,:);
    
    LRP(1).Condition    = [LRP(1).Condition '_' LRP(2).Condition '_LRP'];
    LRP(1).Data         = (LRP(1).Data+LRP(2).Data)/2;
    
    LRP                 = LRP(1);
end

end