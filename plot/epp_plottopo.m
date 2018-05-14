% PURPOSE:  plot topos (averaged across subjects) using EEGLAB's topoplot
%           function.
%
% FORMAT
% ------
% epp_plottopo(study,chanlocs,conditions,timePoints,varargin)
%
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
% chanlocs      - channel locations (as per EEGLAB).
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% timePoints    - a 1-by-times vector of time points to be plotted
%                 (e.g. [0 50 200]).
%               - a 2-by-times matrix with rows indicating the start and
%                 end of the time window to be average across
%                 (e.g., [0 100; 300 350]').
%
% The available parameters are as follows:
%           'plotlabels'    - see topoplot for options (defult: 'on').
%           'maplimits'     - [min max] values for the color scale.
%           'R'             - if true, plot data is saved into a csv file
%                             in the current directory + an R file with
%                             code to plot your topos. (in R you can
%                             continue to format you plot - colors,
%                             annotations, etc.)
%
% See also epp_plotbutterfly, epp_plotgrands, epp_plotTF, epp_plottopoTF, epp_plotchannels
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
14-05-2018  Improvment to exporting plot data
            Silenced topoplot();
01-05-2018  Changed colormap (removed custom option)
15-04-2018  1. Removed smoothing (that wasnet working). Instead added
               support for plotting over explicit time windows.
            2. Removed plotting t-values (was not being used + for
               simplicity) 
            3. Added control of color-map
06-03-2018  Removed samplingRate & baseLine fields from study struct.
05-03-2018  Fix title printing
04-07-2017  Added time-smoothing
16-06-2017  Minor fix for plotting channels in R code
23-04-2017  New function (written in MATLAB R2015a)
%}


function epp_plottopo(study,chanlocs,conditions,timePoints,varargin)


%% Validate Input
p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'chanlocs',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'timePoints',@(x) isnumeric(x) & size(x,1) <= 2 & ndims(x)<=2);
    addParameter(p,'plotlabels','on',@isstr)
    addParameter(p,'maplimits',nan,@isnumeric)
    addParameter(p,'R',false,@islogical)
parse(p, study,chanlocs,conditions,timePoints,varargin{:}); % validate


%% Prepare Data for Plotting


% Get only relevant conditions (in order!)
% ----------------------------------------
cInd  = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study = study(cInd);
clear cInd


% Get 'chanlocs'
% --------------
if isempty(chanlocs) % if missing, load from SET file
    EEG         = pop_loadset;
    chanlocs    = EEG.chanlocs;
end


% Get Time Points
% ---------------
nTimes = size(timePoints,2);
% find time point indecies
timePoints  = sort(timePoints,2); % sort time points
for t = 1:nTimes
    temp_tp1 = dsearchn(study(1).timeLine',timePoints(1,t)); % approximated to existing time points
    temp_tp2 = dsearchn(study(1).timeLine',timePoints(size(timePoints,1),t)); % approximated to existing time points
    tP_ind{t} = temp_tp1:temp_tp2;
end

% Get Values to Plot:
% -------------------
for c = 1:length(study)
    % get time-relevant data only
    for t = 1:nTimes
        % avarage across time window
        temp_data(:,t,:) = mean(study(c).Data(:,tP_ind{t},:),2);
    end
    % average across subjects
    meanData(:,:,c) = mean(temp_data,3); % size = electrodes * time * conditions
    clear temp_data
end

% Get min-max values:
% -------------------
if isnan(p.Results.maplimits)
    maxlim      = max(meanData(:));
    minlim      = min(meanData(:));
    maplimits   = [-max(maxlim, -minlim), max(maxlim, -minlim)];
else
    maplimits = p.Results.maplimits;
end
    


%% Plot
% Each conditions has its own figure.
% Multiple time points are plotted in same figure. 

topo_args = {'style',     'map',... 'both'
            'maplimits',  maplimits,... 'absmax',...
            'electrodes', p.Results.plotlabels, ... 'labels' / 'numpoint'
            'whitebk','on'};

for c = 1:size(meanData,3) % each condition
    fig = figure();
    clf
    hold on;
    set(gca,'Color',get(fig,'Color'));
    
    for t = 1:nTimes % each time point
        subplot(1,nTimes,t);

        % Plot the Data 
        % -------------
        evalc('topoplot(meanData(:,t,c), chanlocs,topo_args{:});');
        
        % Add Color Bar
        % -------------
        if t==1
            cl = colorbar('west','AxisLocation','out');
            cl.Position(1)  = cl.Position(1)*0.4;
            cl.Position(3)  = cl.Position(3)*0.3;
            
            cl.Position(2)  = cl.Position(2)/1.5;
            cl.Position(4)  = cl.Position(4)*1.5;
        end
        
        % Add Title(s)
        % ------------
        if t==1
            suptitle([study(c).Condition]);
        end
        title([num2str(timePoints(:,t)') ' ms'], 'Interpreter', 'none');
        
    end
    colormap(suppMakeColormap('wcbkryw'))
end

%% Save to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    save_data.Condition = {};
    save_data.TimePnt   = {};
    save_data.Channel   = {};
    save_data.Theta     = [];
    save_data.Radius    = [];
    save_data.Amp       = [];

    % Orgenize time points:
    for t = 1:size(timePoints,2)
        temp{t} = num2str(timePoints(:,t)');
    end
    times_str   = temp;
    nTimes      = length(times_str);

    % Orgenize channels
    chanlocs    = struct2table(chanlocs);
    nChans      = size(chanlocs,1);

    for c = 1:length(study) % for each condition
        temp_cond   = repmat({study(c).Condition},	[nChans nTimes])';
        temp_time   = repmat(times_str,             [nChans 1])';
        temp_chan   = repmat(chanlocs.labels,       [1 nTimes])';
        temp_theta  = repmat(chanlocs.theta,        [1 nTimes])';
        temp_radius = repmat(chanlocs.radius,       [1 nTimes])';
        temp_amp    = meanData(:,:,c)';

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.TimePnt   = [save_data.TimePnt;   temp_time(:)];
        save_data.Channel   = [save_data.Channel;   temp_chan(:)];
        save_data.Theta     = [save_data.Theta;     temp_theta(:)];
        save_data.Radius    = [save_data.Radius;    temp_radius(:)];
        save_data.Amp       = [save_data.Amp;       temp_amp(:)];
    end

    T = struct2table(save_data);
    
    % Save to CSV
    % -----------
    fn  = ['topoplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv
    
    
    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\epp_plottopo.R'],'rt');
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