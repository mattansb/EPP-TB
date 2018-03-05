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
%           'all'           - if true, plots all channels across subjects.
%                             (electrodes is ignored.)
%
% See also epp_plotgrands, epp_plotTF, epp_plottopo
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
05-03-2018  Fix title printing
04-03-2018  Added ability to plot Trace plots
18-06-2017  Added ability to plot Jackknifed waves.
18-06-2017  Removed ability to plot time window.
04-04-2017  Added ability to mark subjects on Z axis
25-11-2016  New function (written in MATLAB R2015a)

2DO
- Make option to export to R.
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
parse(p, study, conditions, electrodes, varargin{:}); % validate

%% Get only relevant conditions (in order!)

cInd = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);

study = study(cInd);

if p.Results.all
    [study, nsubs]  = suppMatchSubjects(study,conditions);
    electrodes      = 1:nsubs;
    for c = 1:length(study)
        study(c).Data   = permute(study(c).Data,[3 2 1]);
        study(c).IDs    = table([1:size(study(c).Data,3)]','VariableNames',{'ID'});
    end
elseif p.Results.jackknife
    for c = 1:length(study)
        study(c).Data = suppJackknife('in',study(c).Data);
    end
end

%% Compute aves

for c = 1:length(study)
    study(c).Data   = squeeze(mean(study(c).Data(electrodes,:,:),1));   % Mean across electrodes
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
    color = co(mod(c,length(co)),:);
    
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

end