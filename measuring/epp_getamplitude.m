% PURPOSE:  measure ERP amplitudes.
%
% FORMAT
% ------
% results = epp_getamplitude(measure,study,conditions, electrodes,timeWindow,varargin)
% 
%
% INPUTS
% ------
% measure       - 'peak' / 'mean' / 'integral' / 'area'
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% electrodes    - vector of electrodes to be plotted after averaging (e.g.
%                 [87 85, 92]).
% timeWindow    - two time points ([start end], in ms) within which the
%                 measument will be taken.
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
% See also epp_getlatency
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
06-03-2018  Removed samplingRate & baseLine fields from study struct.
07-02-2017  Added support for data interpolating
19-01-2017  Support for saving data
25-11-2016  New function (written in MATLAB R2015a)

2DO
make point amp measure
%}

function results = epp_getamplitude(measure,study,conditions, electrodes,timeWindow,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
    
    addParameter(p,'jackknife', false, @islogical); % (+fix?) < leave fix for later...
    addParameter(p,'average',false,@islogical);
    addParameter(p,'plot',false,@islogical);
    addParameter(p,'save','no', @ischar);
    addParameter(p,'interpolate',1, @isnumeric);
    
    switch lower(measure)
        % case 'point' %????????????????????????????????????????????
        case 'peak'
            addParameter(p,'direction', 1, @isnumeric);
            addParameter(p,'local', 1, @isnumeric);
        case 'mean'
        case 'integral'
        case 'area'
            addParameter(p,'direction', 0, @isnumeric);
            addParameter(p,'boundary',0,@isnumeric);
        otherwise
            error('no such measure')
    end
parse(p ,study, conditions, electrodes, timeWindow, varargin{:}); % validate

samplingRate = 1000/(study(1).timeLine(2) - study(1).timeLine(1));

%% Orgenize Data Before Measuring

study = suppPrep4stats(study, conditions, electrodes, timeWindow, p.Results.average, p.Results.jackknife, p.Results.interpolate);


%% Get Amplitudes

for c = 1:length(study)
    fprintf('\nCalculating amplitudes for %s (condition %d of %d)..',study(c).Condition, c ,length(study))
    switch lower(measure)
        % case 'point'
        case 'peak'
            [study(c).amplitudes, ~] = suppPeak(study(c).cutData, study(c).cutTime, p.Results.direction, p.Results.local);
        case 'mean'
            study(c).amplitudes = mean(study(c).cutData,2);                     % compute mean
            study(c).amplitudes = permute(study(c).amplitudes,[3 1 2]);         % reorder dimentions
        case 'integral'
            study(c).amplitudes = sum(study(c).cutData,2);                      % compute integral
            study(c).amplitudes = study(c).amplitudes/samplingRate;    % convert to mv/sec
            study(c).amplitudes = permute(study(c).amplitudes,[3 1 2]);         % reorder dimentions
        case 'area'
            study(c).cutData = study(c).cutData-p.Results.boundary;             % offset by boundary

            if p.Results.direction > 0      % only use positive values
                study(c).cutData(study(c).cutData<0) = 0;
            elseif p.Results.direction < 0  % only use negative values
                study(c).cutData(study(c).cutData>0) = 0;
                study(c).cutData = abs(study(c).cutData);
            else                            % use all values
                study(c).cutData = abs(study(c).cutData);
            end

            study(c).amplitudes = sum(study(c).cutData,2);                      % compute area
            study(c).amplitudes = study(c).amplitudes/samplingRate;    % convert to mv/sec
            study(c).amplitudes = permute(study(c).amplitudes,[3 1 2]);         % reorder dimentions
    end
    study(c).measure = study(c).amplitudes;
    fprintf('. Done!')
end

%% Prep for export:

results = suppPrep4export(measure,study,conditions,electrodes,p.Results);

%% Plot Results

if p.Results.plot
    suppPlotResults(study, timeWindow)
end

end