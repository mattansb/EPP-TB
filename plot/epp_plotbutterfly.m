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
%           'all'           - if true, plots selected channels across
%                             subjects (if electrodes is left blank, plots
%                             all channels). 
%
% See also epp_plotgrands, epp_plotTF, epp_plottopo, epp_plottopoTF
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
13-05-2018  Fix to trace plot + added ability to plot trace plots with
            selected channels.
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
    addParameter(p,'all', false, @islogical)
    addParameter(p,'R', false, @islogical)
parse(p, study, conditions, electrodes, varargin{:}); % validate

%% Get only relevant conditions (in order!)

cInd = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);

study = study(cInd);

if p.Results.jackknife % jackknife data
    for c = 1:length(study)
        study(c).Data = suppJackknife('in',study(c).Data);
    end
end

% Select electrodes
if p.Results.all
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

%% Plot
% Find number of minimal subplot dimentions (aprox)
a = ceil(sqrt(length(study)/1.6));
b = ceil(length(study)/a);
if a>b
    b = ceil(sqrt(length(study)/1.6));
    a = ceil(length(study)/b);
end

% Open Figure
fig         = figure();
fig.Color   = [1 1 1];
co          = get(gca,'ColorOrder');
hold on;
clf

for c = 1:length(study) % for each condition
    cond(c) = subplot(a,b,c);                   % open new subplot
    
    if p.Results.minusUp
        set(gca,'YDir','reverse');
        hold on
    end
    
    % Plot the lines
    % --------------
    co_ind = mod(c,length(co));
    if co_ind==0, co_ind = length(co); end
    color = co(co_ind,:);
    
    [nTime, nID]=size(study(c).Data);
    plot3(study(c).timeLine,study(c).Data,repmat(1:nID,[nTime 1]),...
        'Color',color,...
        'LineWidth',0.5);                       % plot "butterflys"
    hold on
    plot(study(c).timeLine,study(c).mean,...
        'Color',[0 0 0],...
        'LineWidth',1.25);                      % plot the mean ERP
    
    set(gca,'ZTickLabel',[study(c).IDs.ID]);

    

    % Limits
    % ------
    xlim([study(c).timeLine([1 end])]); % set x-Axis limit
    
    
    
    % Plot axes
    % ---------
    plot(study(1).timeLine([1,end]), [0 0],...
        'Color', 'k',...
        'LineStyle', '-');                      % plot time line
    plot([0 0],[-2000 2000],...
        'Color', 'k',...
        'LineStyle', '-');                      % plot y-axis (at t=0)
    
    % Titles (and label)
    % ------------------
    title(conditions(c), 'Interpreter', 'none')
    if c == 1 % for first plot only
        xlabel('Time');
        ylabel('\muV')
    end
    
    % change view angle
    view(0,90)
    
end


for c = 1:length(study) % set all plots to same scale
    ylim(cond(c),[min(minA) max(maxA)]); 
end


%% Export to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    nTmPnts = length(study(1).timeLine);    % length of time axis
    nCond   = length(study);                % number of conditions
    nSumID  = sum(cellfun(@(x) size(x,1),{study.IDs}));
    
    % make empty matrices
    ERPamp      = nan(nTmPnts,nSumID);
    ERPcond     = cell(nTmPnts,nSumID);
    ERPID       = cell(nTmPnts,nSumID);

    for c = 1:nCond
        x_start = find(all(isnan(ERPamp),1),1);
        x_end   = x_start + size(study(c).Data,2) - 1;
        
        ERPamp(:,x_start:x_end)     = study(c).Data;
        ERPcond(:,x_start:x_end)    = {study(c).Condition};
        ERPID(:,x_start:x_end)      = table2cell(repmat(study(c).IDs(:,1), 1, nTmPnts))';
    end
    
    ERPtime = repmat(study(1).timeLine',1,nSumID); % get time line
    
    % Reshape all:
    ERPamp  = reshape(ERPamp,1,[]);
    ERPID   = reshape(ERPID,1,[]);
    ERPcond = reshape(ERPcond,1,[]);
    ERPtime = reshape(ERPtime,1,[]);
    
    % Convert to table and write to CSV:
    T   = table(ERPcond', ERPtime', ERPID', ERPamp', 'VariableNames', {'Condition','Time','ID','amp'});
    fn  = ['butterflyplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv
    
    
    % Make R code file
    % ----------------
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