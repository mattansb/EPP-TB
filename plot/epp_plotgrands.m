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
%                             (false my default).
%           'R'             - if true, plot data is saved into a csv file
%                             in the current directory + an R file with
%                             code to plot your ERPs. (in R you can
%                             continue to format you plot - colors,
%                             annotations, etc.)
%
% See also epp_plotbutterfly, epp_plotTF, epp_plottopo, epp_plottopoTF, epp_plotchannels
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
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
parse(p, study, conditions, electrodes, varargin{:}); % validate

errorType = p.Results.errorType;

%% Get only relevant conditions (in order!)

cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study   = study(cInd);

% determine if to compute within
if ~strcmpi(errorType,'NaN')
    warning('When plotting errors, assuming all conditions to be within-subject.')
    nsubs_old = cellfun(@(x) size(x,3), {study.Data});              % how many subjects are in each conditions
    [study_new , nsubs_new] = suppMatchSubjects(study,conditions);  % get data only for subjects that have data in all conditions
    
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
    
    maxA(c) = max(study(c).mean(:));
    minA(c) = min(study(c).mean(:));
end
maxA = max(maxA);
minA = min(minA);

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
    
    % Prep for plotting
    % -----------------
    TT  = [study(1).timeLine, fliplr(study(1).timeLine)];
    for c = 1:length(study)
        switch lower(errorType(1:2))
            case 'sd'
                YY = study(c).sd;
            case 'se'
                YY = study(c).se;
            case 'ci'
                YY = study(c).ci;
        end
        ymin = (study(c).mean - YY)';
        ymax = (study(c).mean + YY)';

        EE(c,:) = [ymin fliplr(ymax)];
    end
    maxA = max(EE(:));
    minA = min(EE(:));
end

%% Plot
fig         = figure();
fig.Color   = [1 1 1];
co          = get(gca,'ColorOrder');
hold on;
clf

for c = 1:length(study)
    co_ind = mod(c,length(co));
    if co_ind==0, co_ind = length(co); end
    color = co(co_ind,:);
    
    % Plot Error(s)
    % =============
    if ~strcmpi(errorType,'NaN')
        errorPlot(c)    = fill(TT, EE(c,:), color,...
                               'EdgeColor', 'none', 'FaceAlpha', 0.25);
    end
    hold on;
    
    % Plot Mean(s)
    % ============
    MM(c) = plot(study(c).timeLine, study(c).mean, 'color', color);
end

% plot Axies
% ==========
plot(study(1).timeLine([1,end]), [0 0],'Color', 'k');   % plot time line
plot([0 0],[minA maxA],'Color', 'k');                   % plot y-axis (at t=0)
ylim([minA maxA]);                                      % set Y limits
xlabel('Time')
ylabel('\muV')
if p.Results.minusUp, set(gca,'YDir','reverse'); end

% Add Legend
legend(MM,{study.Condition}, 'Interpreter', 'none');



%% Export to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    save_data.Condition = {};
    save_data.Time      = [];
    save_data.N         = [];
    save_data.mean      = [];
    save_data.sd        = [];

    for c = 1:length(study) % for each condition
        temp_cond   = repmat({study(c).Condition},  size(study(c).mean));
        temp_N      = repmat(study(c).N,            size(study(c).mean));
        try
            temp_SD = study(c).sd;
        catch
            temp_SD = nan(size(study(c).mean));
        end

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.Time      = [save_data.Time;      study(c).timeLine'];
        save_data.N         = [save_data.N;         temp_N(:)];
        save_data.mean      = [save_data.mean;     study(c).mean(:)];
        save_data.sd        = [save_data.sd;        temp_SD(:)];
    end

    T = struct2table(save_data);
    
    % Save to CSV
    % -----------
    fn  = ['erpplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv
    
    
    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\epp_plotgrands.R'],'rt');
    f = fread(fid,'*char')';
    fclose(fid);
    
    f = strrep(f,'@filename@',[fn '_data.csv']);
    
    fid  = fopen([fn '_code.R'],'wt');
    fprintf(fid,'%s',f);
    fclose(fid);
        
else
    fprintf('NOTE: consiter plotting with ggplot in ''R'', \nNOTE: by setting ''R'' to true.\n')
end

end