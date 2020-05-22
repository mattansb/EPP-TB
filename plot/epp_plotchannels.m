% PURPOSE:  plot data by channel in a grid array / channel positions
%           (averaged across subjects). This is a stripped down version of
%           eeglab's plottopo().
%
% FORMAT
% ------
% epp_plotchannels(study,conditions,electrodes,varargin)
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% electrodes    - vector of electrodes to be plotted after averaging (e.g.
%                 [87 85, 92]). If left blank, all channels will be
%                 plotted.
%
% The available parameters are as follows:
%           'chanlocs'      - eeglab chanlocs struct. If not provided,
%                             channel data will be plotted in a grid array.
%           'R'             - if true, plot data is saved into a csv file
%                             in the current directory + an R file with
%                             code to plot your ERPs. (in R you can
%                             continue to format you plot - colors,
%                             annotations, etc.)
%
% When study is a time-frequency structure, two additional parameters must
% be supplied:
%           'type'          - 'erps' or 'itc'.
%           'freqs'         - matrix of frequencies, with each row
%                             containing a range of frequencies to group
%                             together (1st column is lower limit, 2nd
%                             column is upper limit of each range). e.g.
%                             freqs = [1 3; 4 15; 16 28]; 
%                             A given band is selected as so: [low <= freq <= high]
%
% See also epp_plotbutterfly, epp_plotTF, epp_plottopo, epp_plotgrands
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
21-05-2020  Fix some bugs
            Added support for TF plotting..
14-05-2018  New function (written in MATLAB R2017a)
%}
function epp_plotchannels(study,conditions,electrodes,varargin)

%% Validate input data
p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@isnumeric);
    addParameter(p,'chanlocs', [], @isstruct)    
    addParameter(p,'minusUp', false, @islogical)
    addParameter(p,'R', false, @islogical)    
    if ~isfield(study, 'Data')
        addParameter(p,'freqs', [], @isnumeric)
        addParameter(p,'type', '', @ischar)
    end
parse(p, study, conditions, electrodes, varargin{:}); % validate

chanlocs = p.Results.chanlocs;

%% Get only relevant conditions (in order!)

cInd    = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);
study   = study(cInd);

if ~isfield(study, 'Data')
    study = epp_reshapeTF(study, p.Results.freqs, p.Results.type);
    conditions = {study.Condition};
end


%% prep data
if isempty(electrodes), electrodes = 1:size(study(1).Data,1); end % select channels

for c = 1:length(study)
    plot_data(:,:,c) = mean(study(c).Data(electrodes,:,:),3);
end


% Prep axis data
if isempty(chanlocs)
    chan_labels = arrayfun(@(X) ['Channel ' num2str(X)], electrodes, 'UniformOutput', false);
else
    chanlocs = chanlocs(electrodes);
    chan_w = 0.05;
    chan_h = 0.08; 
    
    % convert theta+radius to x*y
    radianTheta = pi/180*[chanlocs.theta];
    chan_x      = [chanlocs.radius].*sin(radianTheta);
    chan_y      = [chanlocs.radius].*cos(radianTheta);
    
    chan_x = num2cell(chan_x);
    chan_y = num2cell(chan_y);
    [chanlocs.x] = chan_x{:};
    [chanlocs.y] = chan_y{:};
end

%% Plot
if isempty(chanlocs)
    p_channels_grid(plot_data, study(1).timeLine, chan_labels, conditions, ...
        'minusUp', p.Results.minusUp, ...
        'R', p.Results.R);
else
    p_channels_array(plot_data, chanlocs, [chan_w chan_h], study(1).timeLine, conditions,...
        'minusUp', p.Results.minusUp, ...
        'R', p.Results.R);
end

end