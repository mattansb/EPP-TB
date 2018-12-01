% PURPOSE:  measure ERP latencies.
%
% FORMAT
% ------
% results = epp_getLat(measure, study, conditions, electrodes, timeWindow, direction, varargin)
% 
%
% INPUTS
% ------
% measure       - 'peak' / 'relative_criterion' / 'baseline_deviation' /
%                 'absolute_criterion' / 'fractional_area'
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% electrodes    - vector of electrodes to be plotted after averaging (e.g.
%                 [87 85, 92]).
% timeWindow    - two time points ([start end], in ms) within which the
%                 measument will be taken.
% direction     - max / min peak to find (see fractional_area below).
%
% The available parameters are as follows:
%           'jackknife'     - {'off'(default)|'on'|'uncentered'|'unweighted'|'uncentered_unweighted'}
%                             measure using the jackknife technique.
%                             if not 'unweighted', jackknifed means and
%                             values will be weighted by study.IDs.nTrials.
%                             If not 'uncentered', un-jackknifed values
%                             will be recentered around a measurment made
%                             around the grand mean ERP, else will be
%                             recentered around the mean jackknifed values.
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
%           'local'         - if larger than zero, peak is defined as the
%                             largest (/smallest) point which is also:
%                               1. larger (/smaller) than one sample on
%                                  either side.
%                               2. larger (/smaller) then the average of
%                                  the N sample points on either side.
%       for 'relative_criterion'
%           'local'         - same as for peak.
%           'percentage'    - [0-1] Latency is the first point before the
%                             peak that is % of peak amplitude.
%           'first_last'    - {'first' [default] | 'last'} look for the
%                             onset or offset of the component.
%       for 'baseline_deviation'
%           'criterion'     - Latency is the first point to be deviate by
%                             more than X standard deviations calculated on
%                             the baseline. 
%           'baseline'      - length of base line (negative or positive)
%       for 'absolute_criterion'
%           'criterion'     - Latency is the first point to be larger
%                             (smaller) than X mV (positive or negative).
%       for 'fractional_area'
%           'direction'     - If 0 will use rectified area.
%           'percentage'    - [0-1] Latency is the point deviding the area
%                             into %X. Area can be the positive, negative or
%                             rectified area (based on direction [1, -1, 0]).
%           'boundary'      - boundary by which to offset the measurement
%                             of area (in mV).
%
%  ==============================================================
% || When comparing to measurements taken using ERPLAB tool:    ||
% ||    Latency type        Correlation     Ave. error (ms)     ||
% ||    ------------        -----------     ---------------     ||
% ||    local peak          r = 1.00        e = 0.00            ||
% ||    relative_criterion  r = 0.99        e = 0.41            ||
% ||    fractional_area     r = 0.99        e = 1.15            ||
%  ==============================================================
%
% See also epp_getAmp, epp_getTF
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
07-04-2018  New function (written in MATLAB R2017a)
%}
function results = epp_getLat(measure,study,conditions, electrodes,timeWindow,direction,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
    addRequired(p,'direction',@isnumeric);
    
    addParameter(p,'jackknife', 'off', @ischar);
    addParameter(p,'average',false,@islogical);
    addParameter(p,'plot',false,@islogical);
    addParameter(p,'save','no', @ischar);
    addParameter(p,'interpolate',1, @isnumeric);
    
    switch lower(measure)
        case 'peak'
            addParameter(p,'local', 1, @isnumeric);
        case 'relative_criterion' % aka fractional peak
            addParameter(p,'local', 1, @isnumeric);
            addParameter(p,'percentage',0.5,@isnumeric);
            addParameter(p,'first_last','first',@(x) any(strcmp(x,{'first','last'})));
        case 'baseline_deviation' % first pass of X sd (measured in the baseline)
            addParameter(p,'criterion',1,@isnumeric); % required!!
            addParameter(p,'baseline',1,@isnumeric); % required!!
            addParameter(p,'first_last','first',@(x) any(strcmp(x,{'first','last'})));
        case 'absolute_criterion'
            addParameter(p,'criterion',1,@isnumeric); % required!!
            addParameter(p,'first_last','first',@(x) any(strcmp(x,{'first','last'})));
        case 'fractional_area'
            addParameter(p,'percentage',0.5,@isnumeric);
            addParameter(p,'boundary',0,@isnumeric);
        otherwise
            error('no such measure')
    end
parse(p ,study, conditions, electrodes, timeWindow, direction, varargin{:}); % validate

%% Orgenize Data Before Measuring

[study, timeWindow_ind] = suppPrep4statsERP(study, conditions, electrodes, timeWindow,...
    p.Results.average, p.Results.jackknife, p.Results.interpolate);

try
    p.Results.local = p.Results.local*p.Results.interpolate;
end

%% Find latencies

for c = 1:length(study)
    fprintf('\nCalculating latencies for %s (condition %d of %d)..',study(c).Condition, c ,length(study))
    res = nan(size(study(c).Data,2),1);
    for ie = 1:size(study(c).Data,2)
        switch lower(measure)
            case 'peak'
                [~,res(ie)] = m_Peak(study(c).Data(:,ie),timeWindow_ind,direction,...
                    p.Results.local,study(c).timeLine);
            case 'relative_criterion'
                res(ie) = m_latRelative_criterion(study(c).Data(:,ie),timeWindow_ind,direction,...
                    p.Results.local,study(c).timeLine,p.Results.percentage,p.Results.first_last);
            case 'baseline_deviation'
                res(ie) = m_latBaseline_deviation(study(c).Data(:,ie),timeWindow_ind,direction,...
                    p.Results.baseline,study(c).timeLine,p.Results.criterion, p.Results.first_last);
            case 'absolute_criterion'
                res(ie) = m_latAbsolute_criterion(study(c).Data(:,ie),timeWindow_ind,direction,...
                    study(c).timeLine,p.Results.criterion,p.Results.first_last,p.Results.first_last);
            case 'fractional_area'
                res(ie) = m_latFractional_area(study(c).Data(:,ie),timeWindow_ind,direction,...
                    p.Results.boundary,study(c).timeLine,p.Results.percentage);
        end
    end
    study(c).measure = res;
    fprintf('. Done!')
end

%% Prep for export & save(?)

[results, study] = suppPrep4exportERP(measure,study,conditions,electrodes,p.Results);

%% Plot Results

if p.Results.plot
    suppPlotResultsERP(study, timeWindow)
end


end