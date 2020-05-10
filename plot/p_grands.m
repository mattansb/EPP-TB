% PURPOSE:  plot grand average ERP (averaged across subjects and selected
%           electrodes), with or without error terms. 
%
% FORMAT
% ------
% p_grands(grands, errors, timeLine, condLabels, varargin)
% 
%
% INPUTS
% ------
% grands        - time * condition matrix of mean ERPs.
% errors        - time * condition matrix of the ERPs' errors. Drawn
%                 symmetrically around the grands. Can be left empty. 
% timeLine      - Vector of time points for the x-axis, matching 'grands'.
% condLabels    - Cellarry of labels for the conditions, matching 'grands'.
% 
%
%
% The available parameters are as follows:
%       'minusUp'	- if true, plot is flipped so minus is up (false by
%                     default). 
%       'save'      - if true, plot data is saved into a csv file in the
%                     current directory + an R file with code to plot your
%                     ERPs. (in R you can continue to format you plot -
%                     colors, annotations, etc.) 
%       'errorType' - label to give the errors in the saved files.
%
% See also epp_plotgrands
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
10-05-2020  New function (written in MATLAB R2017b)
%}
function p_grands(grands, errors, timeLine, condLabels, varargin)

p = inputParser;
    addRequired(p,'grands',@isnumeric);
    addRequired(p,'errors',@isnumeric);
    addRequired(p,'timeLine',@isnumeric);
    addRequired(p,'condLabels',@iscell);
    addParameter(p,'minusUp',false);
    addParameter(p,'save',false);
    addParameter(p,'errorType','errors',@isstr);
parse(p, grands, errors, timeLine, condLabels, varargin{:}); % validate

minusUp     = p.Results.minusUp;
save        = p.Results.save;
errorType   = p.Results.errorType;

timeemit  = [timeLine, fliplr(timeLine)];


if ~isempty(errors)
    maxA = max(grands(:) + errors(:));
    minA = min(grands(:) - errors(:));
else
    maxA = max(grands(:));
    minA = min(grands(:));    
end
%% Plot
fig         = figure();
fig.Color   = [1 1 1];
co          = get(gca,'ColorOrder');
hold on;
clf

for c = 1:size(grands,2)
    co_ind = mod(c,length(co));
    if co_ind==0, co_ind = length(co); end
    color = co(co_ind,:);
    
    % Plot Error(s)
    % =============
    if ~isempty(errors)
        Emin = grands(:,c) - errors(:,c);
        Emax = grands(:,c) + errors(:,c);
        
        errorPlot(c) = fill(timeemit, [Emin' fliplr(Emax')],...
            color, 'EdgeColor', 'none', 'FaceAlpha', 0.25);
    end
    hold on;
    
    % Plot Mean(s)
    % ============
    MM(c) = plot(timeLine, grands(:,c), 'color', color);
end

% plot Axies
% ==========
plot(timeLine([1,end]), [0 0],'Color', 'k');    % plot time line
plot([0 0],[minA maxA],'Color', 'k');           % plot y-axis (at t=0)
ylim([minA maxA]);                              % set Y limits
xlim(timeLine([1 end]));                        % set X limits
xlabel('Time');
ylabel('\muV');
if minusUp, set(gca,'YDir','reverse'); end

% Add Legend
legend(MM,condLabels, 'Interpreter', 'none');



%% Export to R
if save
    % Save the data in long format
    % ----------------------------
    save_data.Condition     = {};
    save_data.Time          = [];
    save_data.mean          = [];
    save_data.(errorType) 	= [];
    

    for c = 1:size(grands, 2) % for each condition
        temp_cond   = repmat(condLabels(c), size(grands, 1), 1);
        
        if ~isempty(errors)
            temp_er = errors(:,c);
        else
            temp_er = nan(size(grands, 1), 1);
        end

        save_data.Condition     = [save_data.Condition;     temp_cond(:)];
        save_data.Time          = [save_data.Time;          timeLine'];
        save_data.mean          = [save_data.mean;          grands(:,c)];
        save_data.(errorType)   = [save_data.(errorType);   temp_er(:)];
    end

    T = struct2table(save_data);
    
    % Save to CSV
    % -----------
    fn  = ['erpplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv
    
    
    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\p_grands.R'],'rt');
    f = fread(fid,'*char')';
    fclose(fid);
    
    f = strrep(f,'@filename@',[fn '_data.csv']);
    f = strrep(f,'@error@',errorType);
    
    fid  = fopen([fn '_code.R'],'wt');
    fprintf(fid,'%s',f);
    fclose(fid);
        
else
    fprintf('NOTE: consiter plotting with ggplot in ''R'', \nNOTE: by setting ''R'' to true.\n')
end

end