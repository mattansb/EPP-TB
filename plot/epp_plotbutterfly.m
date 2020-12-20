% PURPOSE:  plot butterfly ERPs (averaged across selected electrodes).
%
% FORMAT
% ------
% epp_plotbutterfly(study,conditions, electrodes,varargin)
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% electrodes    - vector of electrodes to be plotted after averaging (e.g.
%                 [87 85, 92]).
%
% The available parameters are as follows:
%           'jackknife'     - if true, will plot jackknifed erps. (defult
%                             to false). 
%           'minusUp'       - if true, plot is flipped so minus is up
%                             (false my default).
%           'trace'         - if true, plots selected channels across
%                             subjects (if electrodes is left blank, plots
%                             all channels). 
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
% See also epp_plotgrands, epp_plotTF, epp_plottopo, epp_plotchannels
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
21-05-2020  Added support for TF plotting.
29-05-2018  Bug fix when jackknifing
14-05-2018  Improvment to exporting plot data
13-05-2018  Fix to trace plot + added ability to plot trace plots with
            selected channels.
            Changed 'all' to 'trace'.
14-04-2018  Fixed error when plotting more than 6 conditions
05-03-2018  Fix title printing
04-03-2018  Added ability to plot Trace plots
18-06-2017  Added ability to plot Jackknifed waves.
18-06-2017  Removed ability to plot time window.
04-04-2017  Added ability to mark subjects on Z axis
25-11-2016  New function (written in MATLAB R2015a)
%}

function epp_plotbutterfly(study,conditions, electrodes,varargin)
%

%% Validate input data
p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@isnumeric);
    addParameter(p,'minusUp', false, @islogical)
    addParameter(p,'jackknife', false, @islogical)
    addParameter(p,'trace', false, @islogical)
    addParameter(p,'all', false, @islogical)
    addParameter(p,'R', false, @islogical)
    if ~isfield(study, 'Data')
        addParameter(p,'freqs', [], @isnumeric)
        addParameter(p,'type', '', @ischar)
    end
parse(p, study, conditions, electrodes, varargin{:}); % validate

%% Get only relevant conditions (in order!)

cInd = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);
study = study(cInd);

if ~isfield(study, 'Data')
    study = epp_reshapeTF(study, p.Results.freqs, p.Results.type);
    conditions = {study.Condition};
end

if p.Results.jackknife % jackknife data
    for c = 1:length(study)
        study(c).Data = f_jackknife('in',study(c).Data,3);
    end
end

% Select electrodes
if p.Results.all || p.Results.trace
    if isempty(electrodes)
        electrodes = 1:size(study(1).Data,1);
    end
    
    for c = 1:length(study)
        study(c).Data   = permute(study(c).Data,[3 2 1]);
        % Mean across subjects
        study(c).Data   = squeeze(mean(study(c).Data(:,:,electrodes),1));
%         study(c).IDs    = table(electrodes','VariableNames',{'ID'});
        study(c).IDs    = table(arrayfun(@num2str,electrodes,'UniformOutput',false)','VariableNames',{'ID'});
        
    end
else
    for c = 1:length(study)
        % Mean across electrodes
        study(c).Data   = squeeze(mean(study(c).Data(electrodes,:,:),1));
    end
end

%% Plot
waves       = {study.Data};
timeLine    = study.timeLine;
condLabels  = {study.Condition};

IDs = {study.IDs};
IDs = cellfun(@(X) X.ID, IDs, 'UniformOutput', false);

p_butterfly(waves, timeLine, condLabels, ...
    'minusUp', p.Results.minusUp,...
    'R',p.Results.R,...
    'IDs', IDs);

end