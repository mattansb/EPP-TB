% PURPOSE:  measure ERP latencies.
%
% FORMAT
% ------
% results = epp_getlatency(measure, study, conditions, electrodes, timeWindow, direction, varargin)
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
%           'local'         - if larger than zero, peak is defined as the
%                             largest (/smallest) point which is also:
%                               1. larger (/smaller) than one sample on
%                                  either side.
%                               2. larger (/smaller) then the average of
%                                  the N sample points on either side.
%       for 'relative_criterion'
%           'local'         - same a for peak.
%           'percentage'    - Latency is the first point before the peak
%                             that is % of peak amplitude.
%       for 'baseline_deviation'
%           'criterion'     - Latency is the first point to be larger
%                             (smaller) than X standard deviations
%                             calculated on the baseline.
%           'baseline'      - length of base line (negative or positive)
%       for 'absolute_criterion'
%           'criterion'     - Latency is the first point to be larger
%                             (smaller) than X mV (positive or negative).
%       for 'fractional_area'
%           'percentage'    - Latency is the point deviding the area into
%                             %X. Area can be positive, negative or
%                             rectified, based on direction (1, -1, 0,
%                             respectively).
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
% See also epp_getamplitude
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
06-03-2018  Removed samplingRate & baseLine fields from study struct.
11-07-2017  Fix bug when measuring local peak with interpolation
07-02-2017  Added support for data interpolating
19-01-2017  Support for saving data
25-11-2016  New function (written in MATLAB R2015a)
%}

function results = epp_getlatency(measure,study,conditions, electrodes,timeWindow,direction,varargin)
% aa = tic;
%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
    addRequired(p,'direction',@isnumeric);
    
    addParameter(p,'jackknife', false, @islogical);
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
        case 'baseline_deviation' % first pass of X sd (measured in the baseline)
            addParameter(p,'criterion',1,@isnumeric); % required!!
            addParameter(p,'baseline',1,@isnumeric); % required!!
        case 'absolute_criterion'
            addParameter(p,'criterion',1,@isnumeric); % required!!
        case 'fractional_area'
            addParameter(p,'area',0,@isnumeric);
            addParameter(p,'percentage',0.5,@isnumeric);
            addParameter(p,'boundary',0,@isnumeric);
        otherwise
            error('no such measure')
    end
parse(p ,study, conditions, electrodes, timeWindow, direction, varargin{:}); % validate

o_direction = direction; % for fractional_area calculations, must save another "backup" of this variable

%% Orgenize Data Before Measuring

study = suppPrep4stats(study, conditions, electrodes, timeWindow, p.Results.average, p.Results.jackknife, p.Results.interpolate);


%% Find latencies
try
    p.Results.local = p.Results.local*p.Results.interpolate;
end

for c = 1:length(study)
    fprintf('\nCalculating latencies for %s (condition %d of %d)..',study(c).Condition, c ,length(study))
    switch lower(measure)
        case 'peak'
            fprintf('.')
            [~, study(c).latencies] = suppPeak(study(c).cutData, study(c).cutTime, direction, p.Results.local); % get amp and latency for suppPeak function
        case 'relative_criterion' % aka fractional peak
            [study(c).amplitudes, study(c).latencies] = suppPeak(study(c).cutData, study(c).cutTime, direction, p.Results.local); % get amp and latency for suppPeak function
            for s = 1:size(study(c).cutData,3)
                fprintf('.')
                for e = 1:size(study(c).cutData,1)
                    try
                        amp                             = study(c).amplitudes(s,e);
                        lat                             = study(c).latencies(s,e);
                        Cr                              = amp*p.Results.percentage;         % compute criterion = %*peak_amp
                        lat                             = find(study(c).cutTime==lat,1);    % convert laency back to index
                        study(c).cutData(e,[lat:end],s) = nan;                              % remove later latencies

                        if direction > 0 % if max amp
                            ind = find(study(c).cutData(e,:,s) < Cr);   % find where amp is smaller than Criterion
                            ind = ind(end);                             % only the LAST one (i.e. first when moving backwards from peak)

                            % Compute linear approximation:
                            A                       = study(c).cutData(e,ind+1,s)-study(c).cutData(e,ind,s);
                            a                       = Cr - study(c).cutData(e,ind,s);
                            B                       = study(c).cutTime(ind+1)-study(c).cutTime(ind);
                            study(c).latencies(s,e) = study(c).cutTime(ind) + B*(a/A); % SAVE LATENCY!
                        elseif direction < 0 % if min amp
                            ind = find(study(c).cutData(e,:,s) > Cr);   % find where amp is larger than Criterion
                            ind = ind(end);                             % only the LAST one (i.e. first when moving backwards from peak)

                            % Compute linear approximation:
                            A                       = study(c).cutData(e,ind,s)-study(c).cutData(e,ind+1,s);
                            a                       = Cr - study(c).cutData(e,ind+1,s);
                            B                       = study(c).cutTime(ind+1)-study(c).cutTime(ind);
                            study(c).latencies(s,e) = study(c).cutTime(ind+1) - B*(a/A); % SAVE LATENCY!
                        end
                    catch
                        study(c).latencies(s,e) = nan;
                    end
                end
            end
        otherwise
            for s = 1:size(study(c).cutData,3)
                fprintf('.')
                for e = 1:size(study(c).cutData,1)
                    try
                        switch lower(measure)
                            case 'baseline_deviation' % first pass of X sd (measured in the baseline)
                                % Find baseline time window
                                if p.Results.baseline < 0
                                    baselineT = study(c).timeLine >= p.Results.baseline & study(c).timeLine < 0;
                                elseif p.Results.baseline > 0
                                    baselineT = study(c).timeLine <= p.Results.baseline & study(c).timeLine > 0;
                                end

                                Cr = std(study(c).Data(e,baselineT,s))*p.Results.criterion; % compute criterion = n*SD
                            case 'absolute_criterion'
                                Cr = p.Results.criterion;
                            case 'fractional_area'
                                % offset amps by boundery
                                study(c).cutData(e,:,s) = study(c).cutData(e,:,s)-p.Results.boundary; 

                                % Rectify area
                                if o_direction == 0 % all area
                                    study(c).cutData(e,:,s)     = abs(study(c).cutData(e,:,s)); % make all amps positive!
                                    direction = 1;
                                elseif o_direction > 0 % only positive area
                                    areaB                       = study(c).cutData(e,:,s) < 0;	% find where amp is negative
                                    study(c).cutData(e,areaB,s) = 0;                            % make neg amps = 0
                                elseif o_direction < 0 % only negative area
                                    areaB                       = study(c).cutData(e,:,s) > 0;  % find where amp is positive
                                    study(c).cutData(e,areaB,s) = 0;                            % make pos amps = 0
                                end

                                % turn amps to cumulative area: 
                                study(c).cutData(e,:,s) = cumsum(study(c).cutData(e,:,s));

                                % compute criterion = %*area
                                Cr = study(c).cutData(e,end,s)*p.Results.percentage;
                        end

                        % find latency!
                        if direction > 0 % if max amp
                            ind = find(study(c).cutData(e,:,s) > Cr,1);   % find first time amp crosses criterion

                            % Compute linear approximation:
                            A                       = study(c).cutData(e,ind+1,s)-study(c).cutData(e,ind,s);
                            a                       = Cr - study(c).cutData(e,ind,s);
                            B                       = study(c).cutTime(ind+1)-study(c).cutTime(ind);
                            study(c).latencies(s,e) = study(c).cutTime(ind) + B*(a/A); % SAVE LATENCY!
                        elseif direction < 0 % if min amp
                            ind = find(study(c).cutData(e,:,s) < Cr,1);   % find first time amp crosses criterion

                            % Compute linear approximation:
                            A                       = study(c).cutData(e,ind,s)-study(c).cutData(e,ind+1,s);
                            a                       = Cr - study(c).cutData(e,ind+1,s);
                            B                       = study(c).cutTime(ind+1)-study(c).cutTime(ind);
                            study(c).latencies(s,e) = study(c).cutTime(ind+1) - B*(a/A); % SAVE LATENCY!
                        end
                    catch
                        study(c).latencies(s,e) = nan;
                    end

                end
            end
    end
    study(c).measure = study(c).latencies;
    fprintf('. Done!')
end

%% Prep for export:

results = suppPrep4export(measure,study,conditions,electrodes,p.Results);

%% Plot Results

if p.Results.plot
    suppPlotResults(study, timeWindow)
end

% toc(aa)
end