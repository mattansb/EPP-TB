% PURPOSE:  plot grand average ERP (averaged across subjects and selected
%           electrodes), with or without error terms. 
%
% FORMAT
% ------
% p_butterfly(waves, timeLine, condLabels, varargin)
% 
%
% INPUTS
% ------
% waves         - Cell list (each cell a condition) of time * waves matrix.
% timeLine      - Vector of time points for the x-axis, matching each cell of 'waves'.
% condLabels    - Cellarry of labels for the conditions, matching 'waves'.
% 
%
%
% The available parameters are as follows:
%       'minusUp'	- if true, plot is flipped so minus is up (false by
%                     default). 
%       'R'         - if true, plot data is saved into a csv file in the
%                     current directory + an R file with code to plot your
%                     ERPs. (in R you can continue to format you plot -
%                     colors, annotations, etc.) 
%       'IDs'       - A cell list (same length as 'waves') of IDs for each
%                     wave. Used when saving the data. 
%
% See also epp_plotbutterfly
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
21-05-2020  Change yaxis to 'value'
10-05-2020  New function (written in MATLAB R2017b)
%}
function p_butterfly(waves, timeLine, condLabels, varargin)

p = inputParser;
    addRequired(p,'waves',@iscell);
    addRequired(p,'timeLine',@isnumeric);
    addRequired(p,'condLabels',@iscell);
    addParameter(p,'minusUp',false);
    addParameter(p,'R',false);
    addParameter(p,'IDs',[]);
parse(p, waves, timeLine, condLabels, varargin{:}); % validate



maxA = max(cellfun(@(X) max(X(:)), waves));
minA = min(cellfun(@(X) min(X(:)), waves));

nConds  = length(waves);
nTime   = length(timeLine);

%% Plot
% Find number of minimal subplot dimentions (aprox)
a = floor(sqrt(nConds));
b = ceil(nConds/a);

% Open Figure
fig         = figure();
fig.Color   = [1 1 1];
co          = get(gca,'ColorOrder');
hold on;
clf

for c = 1:nConds % for each condition
    subplot(a,b,c); % open new subplot
    
    % Select colors:
    co_ind = mod(c,length(co));
    if co_ind==0, co_ind = length(co); end
    color = co(co_ind,:);
    
    % Plot the lines
    % --------------
    nID = size(waves{c},2);
    % plot "butterflys"
    plot3(timeLine,waves{c},repmat(1:nID,[nTime 1]),...
        'Color',color,'LineWidth',0.5);
    hold on
    % plot the mean ERP
    plot(timeLine,mean(waves{c},2),...
        'Color',[0 0 0],'LineWidth',1.25);
    
    % Plot axes
    % ---------
    plot(timeLine([1,end]), [0 0],'Color', 'k');   % plot time line
    plot([0 0],[minA maxA],'Color', 'k');                   % plot y-axis (at t=0)
    
    % Titles, labels, and limits
    % --------------------------
    title(condLabels(c), 'Interpreter', 'none');
    xlim([timeLine([1 end])]); % set x-Axis limit
    ylim([minA maxA]); 
    view(0,90) % change view angle
    if p.Results.minusUp, set(gca,'YDir','reverse'); end
    if c == 1 % for first plot only
        xlabel('Time');
        ylabel('Value');
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

    for c = 1:nConds % for each condition
        temp_cond   = repmat(condLabels(c),  size(waves{c}));
        temp_time   = repmat(timeLine',    [1 size(waves{c},2)]);
        
        if ~isempty(p.Results.IDs)
            temp_ID = repmat(p.Results.IDs{c}',    [size(waves{c},1) 1]);
        else
            IDc     = 1:size(waves{c},2);
            temp_ID = repmat(IDc,    [size(waves{c},1) 1]);
        end

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.Time      = [save_data.Time;      temp_time(:)];
        save_data.ID        = [save_data.ID;        temp_ID(:)];
        save_data.amp       = [save_data.amp;       waves{c}(:)];
    end

    T = struct2table(save_data);

    % Save to CSV
    % -----------
    fn  = ['butterflyplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true);
    
    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\p_butterfly.R'],'rt');
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