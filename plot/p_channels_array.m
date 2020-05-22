% PURPOSE:  plot data by channel in a channel positions array.
%
% FORMAT
% ------
% p_channels_array(plot_data, chanlocs, chan_wh, times, cond_labels, varargin)
% 
%
% INPUTS
% ------
% plot_data     - A channels * time * conditions matrix.
% chanlocs      - A structure with labels (channel label) and x and y
%                 (lower case, marking 2D position!!)
% times         - vector of times.
% cond_labels   - a cell list of the conditions' labels.
%
% The available parameters are as follows:
%       'minusUp'	- if true, plot is flipped so minus is up (false by
%                     default). 
%       'R'         - if true, plot data is saved into a csv file in the
%                     current directory + an R file with code to plot your
%                     ERPs. (in R you can continue to format you plot -
%                     colors, annotations, etc.) 
%
% See also epp_plotchannels
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
22-05-2020  New function (written in MATLAB R2017b)
%}
function p_channels_array(plot_data, chanlocs, times, cond_labels, varargin)

%% Validate input data
p = inputParser;
    addRequired(p,'plot_data',@isnumeric);
    addRequired(p,'chanlocs',@isstruct);
    addRequired(p,'times',@isnumeric);
    addRequired(p,'cond_labels',@iscell);
    
    addParameter(p,'minusUp', false, @islogical)
    addParameter(p,'R', false, @islogical)    
parse(p, plot_data, chanlocs, times, cond_labels, varargin{:}); % validate

%% Params

nChans = size(plot_data, 1);
nTimes = size(plot_data, 2);
nConds = size(plot_data, 3);

[~,ind0] = min(abs(times));

minA = min(plot_data(:));
maxA = max(plot_data(:));

chan_x = [chanlocs.x];
chan_y = [chanlocs.y];

% standardize x,y positions:
chan_x = chan_x - min(chan_x); % subtract min
chan_y = chan_y - min(chan_y);

chan_x = chan_x/max(chan_x); % divide by max
chan_y = chan_y/max(chan_y);

chan_x = 0.1 + chan_x * 0.8;
chan_y = 0.1 + chan_y * 0.8;

% hight and width
% d = dist([chan_x; chan_y]);
% d = tril(d);
% d(d==0) = nan;
% chan_w = min(d(:));
% chan_h = chan_w * 5 / 8;

chan_h = 0.05; 
chan_w = 0.08;

%% Plots

fig         = figure();
fig.Color   = [1 1 1];
hold on;
clf

for ch = 1:nChans
    pos = [chan_x(ch) chan_y(ch) chan_w chan_h];
    ax = axes('Position',pos);

    chan_plot = plot(squeeze(plot_data(ch,:,:)));
    hold on
    text(double(ind0),double(maxA),chanlocs(ch).labels)

    % time & amp axis
    plot([0 nTimes], [0 0],'Color', 'k');
    plot([ind0 ind0],[minA maxA],'Color', 'k');


    ylim([minA maxA]);
    set(gca,'Visible','off')
    if p.Results.minusUp, set(gca,'YDir','reverse'); end
end

legend(chan_plot,cond_labels,...
    'Interpreter', 'none',...
    'Box','off',...
    'Position',[chan_w chan_h*nConds/2 0 0]);

%% Save to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    save_data.Condition = {};
    save_data.Channel   = {};
    save_data.Time      = [];
    save_data.amp       = [];
    save_data.x         = [];
    save_data.y         = [];

    for c = 1:nConds
        temp_cond   = repmat(cond_labels(c),[nTimes nChans]);
        temp_chan   = repmat({chanlocs.labels},[nTimes 1]);
        temp_time   = repmat(times',[1 nChans]);
        temp_amp    = plot_data(:,:,c)';
        temp_x      = repmat(chan_x,[nTimes 1]);
        temp_y      = repmat(chan_y,[nTimes 1]);

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.Channel   = [save_data.Channel;   temp_chan(:)];
        save_data.Time      = [save_data.Time;      temp_time(:)];
        save_data.amp       = [save_data.amp;       temp_amp(:)];
        save_data.x         = [save_data.x;         temp_x(:)];
        save_data.y         = [save_data.y;         temp_y(:)];
    end

    T = struct2table(save_data);

    % Save to CSV
    % -----------
    fn  = ['chanplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv


    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\p_channels_array.R'],'rt');
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