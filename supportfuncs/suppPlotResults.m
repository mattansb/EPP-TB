% This function plots results obtained by the measurement functions - and
% is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
% See also epp_getamplitude & epp_getamplitude

%{
Change log:
-----------
30-01-2017  Changed error bars to denote SE (not 95% CI)
25-11-2016  New function (written in MATLAB R2015a)
%}

function suppPlotResults(study, timeWindow)

%% find minimal number of subplots needed.
% Find number of minimal subplot dimentions (aprox)
a = ceil(sqrt((1+length(study))/1.6));
b = ceil((1+length(study))/a);
if a>b
    b = ceil(sqrt((1+length(study))/1.6));
    a = ceil((1+length(study))/b);
end

%% Plot

% Open Figure
% ===========
fig = figure();
set(fig,'Color',[1 1 1]);
set(fig,'Name','Errors denote SE')
clf
hold on;
set(gca,'Color',get(fig,'Color'));

for c = 1:length(study)
    cond(c) = subplot(a,b,c); % for each subplot:
    hold on;
    
    % Plot Time Window Box /Line
    % ==========================    
    % Get parameters:
    x = timeWindow(1);
    w = timeWindow(2) - timeWindow(1);

    % Plot:
    rectangle('Position',[x -100 w 200],...
        'Curvature',[0 0],...
        'EdgeColor','none',...
        'FaceColor',[0.85 0.85 0.85]);
    
    % Plot Grand ERP
    % ==============
    % Get parameters:
    N(c) = size(study(c).Data,3);
    study(c).Data = squeeze(nanmean(nanmean(study(c).Data(:,:,:),1),3));
    
    % Plot ERP
    plot(study(c).timeLine,study(c).Data,...
        'Color',[0 0 0],...
        'LineWidth',1);

    % Plot axes
    % ---------
    plot(study(1).timeLine([1,end]), [0 0],...
        'Color', 'k',...
        'LineStyle', '-');                      % plot time line
    plot([0 0],[-200 200],...
        'Color', 'k',...
        'LineStyle', '-');                      % plot y-axis (at t=0)
    if c==1
        xlabel('Time');
        ylabel('\muV');
    end
    
    maxA(c) = max(study(c).Data);
    minA(c) = min(study(c).Data);
    
    
    xlim([study(c).timeLine([1,end])]);
    
    title([study(c).Condition]);
end

for c = 1:length(study) % set all plots to same scale
    ylim(cond(c),[min(minA) max(maxA)]);
end



%% Compare Results

% Get parameters
for c = 1:length(study)
    meanValue(c) = nanmean(nanmean(study(c).measure,2),1);
    
    StdE(c) = 0;
    try
%         StdE(c) = nanstd(nanmean(study(c).measure,2),1);
%         Tc = tinv(1-(0.05/2),(N(c)-1));
%         StdE(c) = StdE(c)*(Tc/(sqrt(N(c)-1)));
        StdE(c) = nanstd(nanmean(study(c).measure,2),1)/(sqrt(N(c)-1));
    end
end

% plot
subplot(a,b,a*b)
barChart = bar([meanValue]);
set(barChart,'FaceColor',[0.6 0.6 0.6]);
set(barChart,'EdgeColor',[0 0 0]);
% 
set(gca,'XTickLabel',[{study(:).Condition}]);
hold on;

try
    errorChart = errorbar(1:length(study),meanValue,StdE,'.');    
    set(errorChart,'Color',[0 0 0]);
    set(errorChart,'LineStyle','none');
    set(errorChart,'Marker','none');   
    set(gca,'YLim',[min(meanValue)-max(StdE)*2 max(meanValue)+max(StdE)*2]);
end

title('Result Summary \fontsize{10}(SE)');
end