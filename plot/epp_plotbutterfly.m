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
%
% See also epp_plotgrands, epp_plotTF, epp_plottopo, epp_plottopoTF, epp_plotchannels
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
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
parse(p, study, conditions, electrodes, varargin{:}); % validate

%% Get only relevant conditions (in order!)

cInd = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);

study = study(cInd);

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
        study(c).IDs    = table(electrodes','VariableNames',{'ID'});
    end
else
    for c = 1:length(study)
        % Mean across electrodes
        study(c).Data   = squeeze(mean(study(c).Data(electrodes,:,:),1));
    end
end

%% Compute aves

for c = 1:length(study)
    study(c).mean   = mean(study(c).Data,2);                            % Mean across subjects
    maxA(c)         = max(study(c).Data(:));                            % find max amplidute
    minA(c)         = min(study(c).Data(:));                            % find min amplidute
end
maxA = max(maxA);
minA = min(minA);

%% Plot
% Find number of minimal subplot dimentions (aprox)
a = floor(sqrt(length(study)));
b = ceil(length(study)/a);

nTime = length(study(1).timeLine);

% Open Figure
fig         = figure();
fig.Color   = [1 1 1];
co          = get(gca,'ColorOrder');
hold on;
clf

for c = 1:length(study) % for each condition
    subplot(a,b,c); % open new subplot
    
    % Select colors:
    co_ind = mod(c,length(co));
    if co_ind==0, co_ind = length(co); end
    color = co(co_ind,:);
    
    % Plot the lines
    % --------------
    nID = size(study(c).Data,2);
    % plot "butterflys"
    plot3(study(c).timeLine,study(c).Data,repmat(1:nID,[nTime 1]),...
        'Color',color,'LineWidth',0.5);
    hold on
    % plot the mean ERP
    plot(study(c).timeLine,study(c).mean,...
        'Color',[0 0 0],'LineWidth',1.25);
    
    % Plot axes
    % ---------
    plot(study(1).timeLine([1,end]), [0 0],'Color', 'k');   % plot time line
    plot([0 0],[minA maxA],'Color', 'k');                   % plot y-axis (at t=0)
    
    % Titles, labels, and limits
    % --------------------------
    title(conditions(c), 'Interpreter', 'none');
    xlim([study(c).timeLine([1 end])]); % set x-Axis limit
    ylim([minA maxA]); 
    view(0,90) % change view angle
    if p.Results.minusUp, set(gca,'YDir','reverse'); end
    if c == 1 % for first plot only
        xlabel('Time');
        ylabel('\muV');
    end
end

%% Export to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    save_data.Condition = {};
    save_data.Time      = [];
    save_data.ID        = {};
    save_data.amp       = [];

    for c = 1:length(study) % for each condition
        temp_cond   = repmat({study(c).Condition},  size(study(c).Data));
        temp_time   = repmat(study(c).timeLine',    [1 size(study(c).Data,2)]);
        temp_ID     = repmat(study(c).IDs{:,1}',    [size(study(c).Data,1) 1]);

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.Time      = [save_data.Time;      temp_time(:)];
        save_data.ID        = [save_data.ID;        temp_ID(:)];
        save_data.amp       = [save_data.amp;       study(c).Data(:)];
    end

    T = struct2table(save_data);

    % Save to CSV
    % -----------
    fn  = ['butterflyplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true);
    
    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\epp_plotbutterfly.R'],'rt');
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