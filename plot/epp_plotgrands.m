% PURPOSE:  plot grand average ERP (averaged across subjects and selected
%           electrodes), with or without error terms. 
%
% FORMAT
% ------
% epp_plotgrands(study,conditions,electrodes,varargin)
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
%           'errorType'     - 'SD' / 'SE' / 'CI%%'. type of error to plot
%                             (default: no error). (e.g. 'CI90' will have
%                             an error of CI 95%).
%           'minusUp'       - if true, plot is flipped so minus is up
%                             (false by default).
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
% See also epp_plotbutterfly, epp_plotTF, epp_plottopo, epp_plotchannels
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
21-05-2020  Added support for TF plotting.
10-05-2020  Moved plotting to p_grands
14-05-2018  Improvment to exporting plot data
            Performace improvment
14-04-2018  Fixed error when plotting more than 6 conditions
05-03-2018  Fix title printing
17-05-2017  Fixed bug when exporting to R
23-03-2017  Fixed error when plotting errors with single condition
07-02-2017  Fixed error plotting according to Morey (2008) and r-cookbook
19-01-2017  Added option to export to R + ggplot2 code for plotting
            Removed plotting time window option
26-11-2016  Changed order of plotting error and mean.
25-11-2016  New function (written in MATLAB R2015a)
%}

function epp_plotgrands(study,conditions,electrodes,varargin)


%% Validate input data
p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    addParameter(p,'errorType','NaN',@(x) strcmpi(x,'SE') || strcmpi(x,'SD') || strcmpi(x(1:2),'CI'))
    addParameter(p,'minusUp', false, @islogical)
    addParameter(p,'R', false, @islogical)
    if ~isfield(study, 'Data')
        addParameter(p,'freqs', [], @isnumeric)
        addParameter(p,'type', '', @ischar)
    end
parse(p, study, conditions, electrodes, varargin{:}); % validate

errorType = p.Results.errorType;

%% Get only relevant conditions (in order!)

cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study   = study(cInd);

if ~isfield(study, 'Data')
    study = epp_reshapeTF(study, p.Results.freqs, p.Results.type);
    conditions = {study.Condition};
end

% determine if to compute within
if ~strcmpi(errorType,'NaN')
    warning('When plotting errors, assuming all conditions to be within-subject.')
    nsubs_old = cellfun(@(x) size(x,3), {study.Data});              % how many subjects are in each conditions
    [study_new , nsubs_new] = epp_matchsubjects(study,conditions);  % get data only for subjects that have data in all conditions
    
    if ~all(nsubs_old == nsubs_new) || length(conditions)==1
        warning('Unable to compute within-subject error. Assuming between-subject conditions.')
        within  = false;
    else
        within  = true;
        study   = study_new;
    end
end

%% Compute aves

% Compute Means
% -------------
for c = 1:length(study)
    study(c).Data   = squeeze(mean(study(c).Data(electrodes,:,:),1));    % Mean across electrodes
    study(c).mean   = mean(study(c).Data,2);                    % Mean across subjects
    study(c).N      = size(study(c).Data,2);                    % number of subjects
end

if ~strcmpi(errorType,'NaN')
    % Compute SD
    % ----------
    if within
        sd_correction = (length(study)/(length(study)-1));  % correction acording to Morey (2008)
        
        % get all data in the same matrix [1, time, subject, condition]
        for c = 1:length(study)
            all_data(:,:,c) = study(c).Data;
        end
        sub_mean    = mean(all_data,3); % subjects' average across conditions (for each time point)
        time_mean   = mean(sub_mean,2); % average from each time point

        % re-size to size of 'all_data'
        sub_mean    = repmat(sub_mean,[1 1 length(study)]);
        time_mean   = repmat(time_mean,[1 study(1).N length(study)]);

        % Compute normed data
        norm_data = all_data - sub_mean + time_mean;

        for c = 1:length(study)
            study(c).sd = std(norm_data(:,:,c),0,2) * sd_correction; % Stds across subject
        end
    else
        for c = 1:length(study)
            study(c).sd = std(study(c).Data,0,2); % Stds across subject
        end
    end
    
    % Compute SE and CI
    % -----------------
    for c = 1:length(study)
        study(c).se = study(c).sd / sqrt(study(c).N);
        if strcmpi(errorType(1:2),'ci')
            % Get alpha and T value
            alpha       = 1 - str2double(errorType(3:end))/100;
            Tc          = tinv((1-alpha/2), study(c).N - 1);    % find critical T
            study(c).ci = study(c).se * Tc;                     % compute CI
        end
    end
    
    % Prep Errors
    % -----------
    for c = 1:length(study)
        switch lower(errorType(1:2))
            case 'sd'
                EE(c,:) = study(c).sd;
            case 'se'
                EE(c,:) = study(c).se;
            case 'ci'
                EE(c,:) = study(c).ci;
        end
    end
end

%% Plot

grands      = cat(2,study.mean);
timeLine    = study(1).timeLine;
condLabels  = {study.Condition};

if ~strcmpi(errorType,'NaN')
    errors = EE';
else
    errors = [];
end

p_grands(grands, errors, timeLine, condLabels, ...
    'minusUp', p.Results.minusUp, ...
    'save', p.Results.R,...
    'errorType', p.Results.errorType)

end