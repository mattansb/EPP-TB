% PURPOSE:  plot data by channel in a grid.
%
% FORMAT
% ------
% p_channels_grid(plot_data, times, chan_labels, cond_labels, varargin)
% 
%
% INPUTS
% ------
% plot_data     - A channels * time * conditions matrix.
% times         - vector of times.
% chan_labels   - a cell list of the channels' labels.
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
function p_channels_grid(plot_data, times, chan_labels, cond_labels, varargin)

%% Validate input data
p = inputParser;
    addRequired(p,'plot_data',@isnumeric);
    addRequired(p,'times',@isnumeric);
    addRequired(p,'chan_labels',@iscell);
    addRequired(p,'cond_labels',@iscell);
    
    addParameter(p,'minusUp', false, @islogical)
    addParameter(p,'R', false, @islogical)    
parse(p, plot_data, times, chan_labels, cond_labels, varargin{:}); % validate

%% Params

nChans = size(plot_data, 1);
nTimes = size(plot_data, 2);
nConds = size(plot_data, 3);
a = floor(sqrt(nChans));
b = ceil(nChans/a);

minA = min(plot_data(:));
maxA = max(plot_data(:));

%% Plots

fig         = figure();
fig.Color   = [1 1 1];
hold on;
clf

for ch = 1:nChans
    
    subplot(a,b,ch); % for each subplot:

    chan_plot = plot(times,squeeze(plot_data(ch,:,:)));

    hold on
    plot(times([1 end]), [0 0],'Color', 'k');
    plot([0 0],[minA maxA],'Color', 'k');

    title(chan_labels{ch});
    ylim([minA maxA]);
    xlim(times([1 end]))
    if p.Results.minusUp, set(gca,'YDir','reverse'); end
end

legend(chan_plot,cond_labels,...
    'Interpreter', 'none',...
    'Box','off',...
    'Orientation','horizontal',...
    'Position',[0 0 1 0.1]);

%% Save to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    save_data.Condition = {};
    save_data.Channel   = {};
    save_data.Time      = [];
    save_data.amp       = [];

    for c = 1:nConds
        temp_cond   = repmat(cond_labels(c),[nTimes nChans]);
        temp_chan   = repmat(chan_labels,[nTimes 1]);
        temp_time   = repmat(times',[1 nChans]);
        temp_amp    = plot_data(:,:,c)';

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.Channel   = [save_data.Channel;   temp_chan(:)];
        save_data.Time      = [save_data.Time;      temp_time(:)];
        save_data.amp       = [save_data.amp;       temp_amp(:)];
    end

    T = struct2table(save_data);

    % Save to CSV
    % -----------
    fn  = ['chanplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv


    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\p_channels_grid.R'],'rt');
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