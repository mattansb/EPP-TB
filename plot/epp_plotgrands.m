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
% See also epp_plotbutterfly, epp_plottopo, epp_plotTF
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
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

% if length(conditions)==1
%     errorType = 'NaN';
% else
    errorType = p.Results.errorType;
% end

%% Get only relevant conditions (in order!)

% cInd    = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);
cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study   = study(cInd);

% If plotting errors, match subjects
% if ~strcmpi(errorType,'NaN')
within = false;
if ~strcmpi(errorType,'NaN') && length(conditions)~=1
    warning('When plotting errors, assuming all conditions to be within-subject.')
    
    nS_old      = cellfun(@(x) size(x,3), {study.Data});        % how many subjects were in each conditions
    
    study_new   = suppMatchSubjects(study,conditions);          % get data only for subjects that have data in all conditions
    
    nS_new      = cellfun(@(x) size(x,3), {study_new.Data});    % how many subjects are now in each conditions
    nSdiff      = nS_old - nS_new;                              % get difference in number of subjects
    
    if any(nS_new == 0)
        warning('Unable to compute within-subject error. Assuming between-subject conditions.')
        within = false;
    else
        if any(nSdiff ~= 0)
            warning('Some subjects were removed for not having data in all conditions')
        end
        within  = true;
        study   = study_new;
    end
end

%% Compute aves

for c = 1:length(study)
    study(c).Data   = mean(study(c).Data(electrodes,:,:),1);    % Mean across electrodes
    study(c).mean   = mean(study(c).Data,3);                    % Mean across subjects
    study(c).N      = size(study(c).Data,3);                    % number of subjects
end

%% Compute Errors
if ~strcmpi(errorType,'NaN')
    if within
        % get all data in the same matrix [1, time, subject, condition]
        for c = 1:length(study)
            all_data(1,:,:,c) = study(c).Data;
        end
        sub_mean    = mean(all_data,4); % subjects' average across conditions (for each time point)
        time_mean   = mean(sub_mean,3); % average from each time point

        % re-size to size of 'all_data'
        sub_mean    = repmat(sub_mean,1,1,1,length(study));
        time_mean   = repmat(time_mean,1,1,study(1).N,length(study));

        % Compute normed data
        norm_data = all_data - sub_mean + time_mean;

        for c = 1:length(study)
            study(c).Data_norm  = norm_data(:,:,:,c);
            study(c).sd         = std(study(c).Data_norm,0,3);                      % Stds across subject
            study(c).sd         = study(c).sd * (length(study)/(length(study)-1));  % correction acording to Morey (2008)
        end
    else
        for c = 1:length(study)
            study(c).sd = std(study(c).Data,0,3); % Stds across subject
        end
    end
    
    for c = 1:length(study)
        study(c).se = study(c).sd/sqrt(study(c).N);
        if strcmpi(errorType(1:2),'ci')
            % Get alpha
            alpha = str2num(errorType(3:end))/100;

            for c = 1:length(study)
                Tc          = tinv((0.5 + alpha/2), study(c).N - 1);    % find critical T
                study(c).ci = study(c).se * Tc;                         % compute CI
            end
        end
    end
end

%% Plot
fig         = figure();
fig.Color   = [1 1 1];
co          = get(gca,'ColorOrder');
hold on;
clf

if p.Results.minusUp
    set(gca,'YDir','reverse');
end

maxAmp = [];
minAmp = [];

for c = 1:length(study)
    color = co(mod(c,length(co)),:);
    
    % Plot Error(s)
    % =============
    if ~strcmpi(errorType,'NaN')
        XX  = [study(c).timeLine, fliplr(study(c).timeLine)];
        MM  = [study(c).mean, fliplr(study(c).mean)];
        
        switch lower(errorType(1:2))
            case 'sd'
                YY = study(c).sd;
            case 'se'
                YY = study(c).se;
            case 'ci'
                YY = study(c).ci;
        end
        
        YY  = [YY, fliplr(YY*(-1))];
        YY = MM + YY;
        
        errorPlot(c)    = fill(XX, YY, color,...
                               'EdgeColor', 'none', 'FaceAlpha', 0.25);
        maxAmp(end+1)   = max(YY);
        minAmp(end+1)   = min(YY);
    end
    hold on;
    
    % Plot Mean(s)
    % ============
    meanPlot(c)     = plot(study(c).timeLine, study(c).mean, 'color', color);
    maxAmp(end+1)   = max(study(c).mean);
    minAmp(end+1)   = min(study(c).mean);
    
end

% plot x-Axes
% ===========
plot(study(1).timeLine([1,end]), [0 0],...
    'Color', 'k',...
    'LineStyle', '-');                  % plot time line
xlabel('Time')                          % label Y-axis

% plot y-Axis
% ===========
plot([0 0],[-2000 2000],...
    'Color', 'k', 'LineStyle', '-');	% plot y-axis (at t=0)
ylabel('\muV')                          % label Y-axis
ylim(1.1*[min(minAmp) max(maxAmp)]);    % set Y limits

% Add Legend
legend([meanPlot],{study.Condition}, 'Interpreter', 'none');



%% Export to R
if p.Results.R
    
    % Save the data in long format
    % ----------------------------
    nTmPnts = length(study(1).timeLine);    % length of time axis
    nCond   = length(study);                % number of conditions
    
    % make empty matrices
    ERPmean     = nan(nTmPnts,nCond);
    ERPsd       = nan(nTmPnts,nCond);
    ERPn        = nan(nTmPnts,nCond);
    ERPcond     = cell(nTmPnts,nCond);
    
    for c = 1:nCond
        ERPmean(:,c)    = study(c).mean';                           % get Mean
        try ERPsd(:,c)  = study(c).sd'; end                         % get SD
        ERPn(:,c)       = repmat(study(c).N, nTmPnts,1);            % get N
        ERPcond(:,c)    = repmat({study(c).Condition},nTmPnts,1);   % get condition name
    end
    
    ERPtime = repmat(study(1).timeLine',1,nCond); % get time line
    
    % Reshape all:
    ERPmean = reshape(ERPmean,1,[]);
    ERPsd   = reshape(ERPsd,1,[]);
    ERPn    = reshape(ERPn,1,[]);
    ERPcond = reshape(ERPcond,1,[]);
    ERPtime = reshape(ERPtime,1,[]);
    
    % Convert to table and write to CSV:
    T   = table(ERPcond', ERPtime', ERPn', ERPmean', ERPsd', 'VariableNames', {'Condition','Time','N','mean','sd'});
    fn  = ['erpplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv
    
    
    % Make R code file
    % ----------------
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