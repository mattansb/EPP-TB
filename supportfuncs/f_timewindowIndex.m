% PURPOSE:  get time window index for use with measurment functions
%
% FORMAT
% ------
% t_ind = f_timewindowIndex(timeWindow, varargin)
%
% EXAMPLES
% --------
% t_ind = f_timewindowIndex([100 150], 'times', -200:4:1000);
% t_ind = f_timewindowIndex([100 150], 'sRate', 250, 'timeRange', [-200 1000]);
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
%{
Change log:
-----------
07-04-2018  New function (written in MATLAB R2017a)
%}
function t_ind = f_timewindowIndex(timeWindow, varargin)

%% Validate

p = inputParser;
    addRequired(p,'timeWindow',@(x) length(x)<=2 && isnumeric(x));
    addParameter(p,'times', nan, @(x) length(x)>1 && isnumeric(x));
    addParameter(p,'sRate', nan, @(x) length(x)==1 && isnumeric(x));
    addParameter(p,'timeRange', nan, @(x) length(x)==2 && isnumeric(x));
parse(p ,timeWindow, varargin{:}); % validate

times       = p.Results.times;
sRate       = p.Results.sRate;
timeRange   = p.Results.timeRange;

%% Set times
if isnan(times)
    if isnan(sRate) || isnan(timeRange(1))
        error('Must supply times or sRate and timeRange.')
    end
    
    times = timeRange(1):(1000/sRate):timeRange(2);
end


%% Get index

t_ind = dsearchn(times',timeWindow'); % closest point to start/end of defined time window
if length(timeWindow)==2
    t_ind = t_ind(1):t_ind(2);
end


end