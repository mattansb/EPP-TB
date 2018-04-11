% PURPOSE:  measure ERP amplitudes.
%
% FORMAT
% ------
% results = epp_getAmp(measure,study,conditions, electrodes,timeWindow,varargin)
% 
%
% INPUTS
% ------
% measure       - 'point' / 'peak' / 'mean' / 'integral' / 'area'
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% electrodes    - vector of electrodes to be plotted after averaging (e.g.
%                 [87 85, 92]).
% timeWindow    - two time points ([start end], in ms) within which the
%                 measument will be taken. If measure is 'point' only one
%                 time point.
%
% The available parameters are as follows:
%           'jackknife'     - measure using the jackknife technique?
%                             (default: false)
%           'average'       - average across electrodes before measuring?
%                             (false my default). 
%           'plot'          - plot results when done? (default: false)
%           'save'          - 'long' / 'wide'; will save the results in the
%                             current directory in the specified format.
%           'interpolate'   - rate by with sampling rate should be
%                             scaled. e.g. is original is 250hz, and
%                             'interpolate' is 2, resulting sampling rate
%                             will be 500hz. Interpolation is done using
%                             the 'spline' method.
%
%       for 'peak'
%           'direction'     - Area can be positive, negative or rectified,
%                             based on direction (1, -1, 0, respectively).
%           'local'         - if larger than zero, peak is defined as the
%                             largest (/smallest) point which is also:
%                               1. larger (/smaller) than one sample on
%                                  either side.
%                               2. larger (/smaller) then the average of
%                                  the N sample points on either side.
%
%       for 'area'
%           'direction'     - Area can be positive, negative or rectified,
%                             based on direction (1, -1, 0, respectively).
%           'boundary'      
%
%  ==============================================================
% || When comparing to measurements taken using ERPLAB tool:    ||
% ||    Latency type    Correlation     Ave. error              ||
% ||    ------------    -----------     ----------              ||
% ||    local peak      r = 0.99        e = 0.00                ||
% ||    mean            r = 0.99        e = 0.00                ||
% ||    integral        r = 0.99        e = 0.00                ||
% ||    area            r = 0.99        e = 0.00                ||
%  ==============================================================
%
% See also epp_getLat, epp_getTF
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
07-04-2018  New function (written in MATLAB R2017a)

%}

function results = epp_getAmp(measure,study,conditions, electrodes,timeWindow,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    
    addParameter(p,'jackknife', false, @islogical); % (+fix?) < leave fix for later...
    addParameter(p,'average',false,@islogical);
    addParameter(p,'plot',false,@islogical);
    addParameter(p,'save','no', @ischar);
    addParameter(p,'interpolate',1, @isnumeric);
    
    switch lower(measure)
        case 'point'
            addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==1);
        case 'peak'
            addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
            addParameter(p,'direction', 1, @isnumeric);
            addParameter(p,'local', 1, @isnumeric);
        case 'mean'
            addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
        case 'integral'
            addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
        case 'area'
            addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
            addParameter(p,'direction', 0, @isnumeric);
            addParameter(p,'boundary',0,@isnumeric);
        otherwise
            error('no such measure')
    end
parse(p ,study, conditions, electrodes, timeWindow, varargin{:}); % validate

samplingRate = 1000/(study(1).timeLine(2) - study(1).timeLine(1));

%% Prep data

[study, timeWindow_ind] = suppPrep4statsERP(study, conditions, electrodes, timeWindow,...
    p.Results.average, p.Results.jackknife, p.Results.interpolate);

try
    p.Results.local = p.Results.local*p.Results.interpolate;
end

%% Get Amplitudes

for c = 1:length(study)
    fprintf('\nCalculating amplitudes for %s (condition %d of %d)..',study(c).Condition, c ,length(study))
    res = nan(size(study(c).Data,2),1);
    for ie = 1:size(study(c).Data,2)
        switch lower(measure)
            case 'point'
                res(ie) = m_ampPoint(study(c).Data(:,ie),timeWindow_ind);
            case 'peak'
                [res(ie),~] = m_Peak(study(c).Data(:,ie),timeWindow_ind,p.Results.direction,p.Results.local,study(c).timeLine);
            case 'mean'
                res(ie) = m_ampMean(study(c).Data(:,ie),timeWindow_ind);
            case 'integral'
                res(ie) = m_ampIntegral(study(c).Data(:,ie),timeWindow_ind,samplingRate);
            case 'area'
                res(ie) = m_ampArea(study(c).Data(:,ie),timeWindow_ind,p.Results.direction,samplingRate,p.Results.boundary);
        end
    end
    study(c).measure = res;
    clear res
    fprintf('. Done!')
end

%% Prep for export & save(?)

[results, study] = suppPrep4exportERP(measure,study,conditions,electrodes,p.Results);

%% Plot Results

if p.Results.plot
    suppPlotResults(study, timeWindow)
end


end