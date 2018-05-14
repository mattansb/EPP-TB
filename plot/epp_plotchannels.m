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
% See also epp_plotbutterfly, epp_plotTF, epp_plottopo, epp_plottopoTF, epp_plotgrands
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
14-05-2018  New function (written in MATLAB R2017a)
%}
function epp_plotchannels(study,conditions,electrodes,varargin)

%% Validate input data
p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@isnumeric);
    addParameter(p,'chanlocs', false, @isstruct)    
    addParameter(p,'R', false, @islogical)    
parse(p, study, conditions, electrodes, varargin{:}); % validate


%% Tidy

cInd = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);

study = study(cInd);

% select channels
if isempty(electrodes)
    electrodes = 1:size(study(1).Data,1);
end

%% prep data
for c = 1:length(study)
    plot_data(:,:,c) = mean(study(c).Data,3);
end
minA = min(plot_data(:));
maxA = max(plot_data(:));

% Prep axis data
chanlocs = p.Results.chanlocs;
if isempty(chanlocs)
    nChans  = length(electrodes);
    a       = floor(sqrt(nChans));
    b       = ceil(nChans/a);
    
    legend_arg = {'Orientation','horizontal',...
        'Position',[0 0 1 0.1]};
else
    chanlocs = chanlocs(electrodes);
    chan_w = 0.05;
    chan_h = 0.08; 
    
    % convert theta+radius to x*y
    radianTheta = pi/180*[chanlocs.theta];
    chan_x      = [chanlocs.radius].*sin(radianTheta);
    chan_y      = [chanlocs.radius].*cos(radianTheta);
    
    % standardize x,y positions:
    chan_x = chan_x - min(chan_x) + chan_w; % subtract min
    chan_y = chan_y - min(chan_y) + chan_h;

    chan_x = chan_x/(max(chan_x)+chan_w); % divide by max
    chan_y = chan_y/(max(chan_y)+chan_h);
    
    % define time axis
    [~,ind0] = min(abs(study(1).timeLine));
    nTimes = length(study(1).timeLine);
    
    legend_arg = {'Position',[chan_w chan_h*length(conditions)/2 0 0]};
end

%% Plot
fig         = figure();
fig.Color   = [1 1 1];
hold on;
clf

for ch = electrodes
    if isempty(chanlocs)
        subplot(a,b,ch); % for each subplot:
        
        chan_plot = plot(study(1).timeLine,squeeze(plot_data(ch,:,:)));
        
        hold on
        plot(study(1).timeLine([1 end]), [0 0],'Color', 'k');
        plot([0 0],[minA maxA],'Color', 'k');

        title(chanlocs(ch).labels);
        ylim([minA maxA]);
        xlim(study(1).timeLine([1 end]))
    else
        pos = [chan_x(ch) chan_y(ch)] - [chan_w chan_h]/2;
        pos = [pos chan_w chan_h];
        ax = axes('Position',pos);

        chan_plot = plot(squeeze(plot_data(ch,:,:)));
        hold on
        text(ind0,maxA,chanlocs(ch).labels)

        % time & amp axis
        plot([0 nTimes], [0 0],'Color', 'k');
        plot([ind0 ind0],[minA maxA],'Color', 'k');


        ylim([minA maxA]);
        set(gca,'Visible','off')
    end
end

legend(chan_plot,conditions,...
    'Interpreter', 'none',...
    'Box','off',...
    legend_arg{:});

%% Save to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    save_data.Condition = {};
    save_data.Channel   = {};
    save_data.Time      = [];
    save_data.amp       = [];
    if ~isempty(chanlocs)
        save_data.x = [];
        save_data.y = [];
    end

    for c = 1:length(conditions)
        temp_cond   = repmat(conditions(c),[length(study(c).timeLine) length(chanlocs)]);
        temp_chan   = repmat({chanlocs.labels},[length(study(c).timeLine) 1]);
        temp_time   = repmat(study(c).timeLine',[1 length(chanlocs)]);
        temp_amp    = plot_data(electrodes,:,c)';

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.Channel   = [save_data.Channel;   temp_chan(:)];
        save_data.Time      = [save_data.Time;      temp_time(:)];
        save_data.amp       = [save_data.amp;       temp_amp(:)];

        try
            temp_x = repmat(chan_x,[length(study(c).timeLine) 1]);
            temp_y = repmat(chan_y,[length(study(c).timeLine) 1]);
            save_data.x = [save_data.x; temp_x(:)];
            save_data.y = [save_data.y; temp_y(:)];
        end
    end

    T = struct2table(save_data);

    % Save to CSV
    % -----------
    fn  = ['chanplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv


    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\epp_plotchannels.R'],'rt');
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