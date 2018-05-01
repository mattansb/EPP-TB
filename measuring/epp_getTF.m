% PURPOSE:  measure ERP amplitudes.
%
% FORMAT
% ------
% results = epp_getTF(measure,study,conditions, electrodes,timeWindow,freqs,varargin)
% 
%
% INPUTS
% ------
% measure       - 'ersp' / 'itc'
% study         - structure.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% electrodes    - vector of electrodes to be plotted after averaging (e.g.
%                 [87 85, 92]).
% timeWindow    - two time points ([start end], in ms) within which the
%                 measument will be taken.
% freqs         - matrix of frequencies, with each row containing a range
%                 of frequencies to group together (1st column is lower
%                 limit, 2nd column is upper limit of each range). e.g.
%                 freqs = [1 3; 4 15; 16 28];
%
% The available parameters are as follows:
%           'average'       - average across electrodes before measuring?
%                             (false my default). 
%           'jackknife'     - measure using the jackknife technique?
%                             (default: false)
%           'save'          - 'long' / 'wide'; will save the results in the
%                             current directory in the specified format.
%
% See also epp_getAmp, epp_getLat
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
%{
Change log:
-----------
18-04-2018  Fix bug when sending to suppPrep4statsTF
14-04-2018  Support for jackknife
07-04-2018  Rewrite of function.
28-02-2018  ITC (abs) is now computed in function, to allow for the
            combination of conditions.
17-08-2017  New function (written in MATLAB R2015a)

2DO
===
p.Results.plot??
%}
function results = epp_getTF(measure,study,conditions, electrodes,timeWindow,freqs,varargin)

%% Validate

p = inputParser;
    addRequired(p,'measure',@(x) any(strcmpi(x,{'ersp','itc'})));
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
    addRequired(p,'freqs',@(x) isnumeric(x) && size(x,2)==2);
    
    addParameter(p,'average',false,@islogical);
    addParameter(p,'jackknife',false,@islogical);
    addParameter(p,'plot',false,@islogical);
    addParameter(p,'save','no', @ischar);
parse(p, measure ,study, conditions, electrodes, timeWindow, freqs, varargin{:}); % validate


%% Orgenize Data Before Measuring

[study, timeWindow_ind, freqs_name] = suppPrep4statsTF(study, conditions, electrodes,...
    timeWindow, freqs,...
    p.Results.average, p.Results.jackknife);

%% Get Measure 

for c = 1:length(study)
    fprintf('\nCalculating for %s (condition %d of %d)..',study(c).Condition, c ,length(study))
    res = nan(size(study(c).ersp,2),1);
    for ie = 1:size(study(c).ersp,2)
        switch measure
            case 'ersp'
                res(ie) = m_tfERSP(study(c).ersp(:,ie),timeWindow_ind);
            case 'itc'
                res(ie) = m_tfITC(study(c).itc(:,ie),timeWindow_ind);
        end
    end
    study(c).measure = res;
    clear res
    fprintf('. Done!')
end

%% Prep for export & save(?)

results = suppPrep4exportTF(measure,study,conditions,electrodes,freqs_name,p.Results);

end