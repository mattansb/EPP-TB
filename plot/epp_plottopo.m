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
% timePoints    - vector of time points to be plotted after averaging (e.g.
%                 [0 50 200]).
%
% The available parameters are as follows:
%           't'             - 'paired' / 'twosample'. Plot t values between two
%                             conditions. 
%           'plotlabels'    - see topoplot for options (defult: 'on').
%           'maplimits'     - [min max] values for the color scale.
%           'R'             - if true, plot data is saved into a csv file
%                             in the current directory + an R file with
%                             code to plot your topos. (in R you can
%                             continue to format you plot - colors,
%                             annotations, etc.)
%           'smooth'        - {num,'ms'} | {num,'sample'} avarage on each
%                             side of each time point. 
%
% See also epp_plotbutterfly, epp_plotgrands, epp_plotTF, epp_plottopoTF
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
06-03-2018  Removed samplingRate & baseLine fields from study struct.
05-03-2018  Fix title printing
04-07-2017  Added time-smoothing
16-06-2017  Minor fix for plotting channels in R code
23-04-2017  New function (written in MATLAB R2015a)

BUGS?
----?
bug when saving t data?
%}


function epp_plottopo(study,chanlocs,conditions,timePoints,varargin)


%% Validate Input
p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'chanlocs',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'timePoints',@isvector);
    addParameter(p,'t','off',@(x) any(strcmpi(x,{'paired','twosample'})) && (length(conditions)==2 || length(conditions)==1))
    addParameter(p,'plotlabels','on',@isstr)
    addParameter(p,'maplimits',nan,@isnumeric)
    addParameter(p,'R',false,@islogical)
    addParameter(p,'smooth',{0},@(x) isnumeric(x{1}) && any(strcmpi(x{2},{'ms','sample','samples'})))
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
% Smoothing factor
nTimes = size(timePoints,2);
if p.Results.smooth{1}~=0
    switch p.Results.smooth{2}
        case 'ms'           
            smooth = p.Results.smooth{1};
        otherwise
            smooth = p.Results.smooth{1}*(study(1).timeLine(2) - study(1).timeLine(1));
    end
    timePoints(2,:) = timePoints(1,:)-smooth;
    timePoints(3,:) = timePoints(1,:)+smooth;
    timePoints      = timePoints(2:3,:);
    clear smooth
end

% find time point indecies
timePoints  = sort(timePoints,2); % sort time points
for t = 1:size(timePoints,1)
    tP_ind(t,:) = dsearchn(study(1).timeLine',timePoints(t,:)')'; % approximated to existing time points
end

% get time-relevant data only
for c = 1:length(study)
    if size(timePoints,1)==2
        for t = 1:nTimes
            temp_data(:,t,:) = mean(study(c).Data(:,tP_ind(1,t):tP_ind(2,t),:),2); % avarage across time window
        end
    else
        temp_data(:,:,:) = study(c).Data(:,tP_ind,:); 
    end
    study(c).Data = temp_data; % size = electrodes * timeP * subs
    clear temp_data
end


% Get Values to Plot:
% -------------------
if ~strcmpi(p.Results.t,'off')
    if strcmpi(p.Results.t,'paired')
        study = suppMatchSubjects(study,conditions);
        [~,~,~,stats] = ttest(study(1).Data,study(2).Data,'dim',3); % paired
    else
        [~,~,~,stats] = ttest2(study(1).Data,study(2).Data,'dim',3); % not paired
    end

    meanData = stats.tstat;
else
    for c = 1:length(study)
        meanData(:,:,c) = mean(study(c).Data,3); % size = electrodes * time * conditions
    end
end


% Get min-max values:
% -------------------
if isnan(p.Results.maplimits)
    maxlim      = max(meanData(:));
    minlim      = min(meanData(:));
    maplimits   = [-max(maxlim, -minlim), max(maxlim, -minlim)];
    clear maxlim minlim
else
    maplimits = p.Results.maplimits;
end
    


%% Plot
% Each conditions has its own figure.
% Multiple time points are plotted in same figure. 

for c = 1:size(meanData,3) % each condition
    fig = figure();
    clf
    hold on;
    set(gca,'Color',get(fig,'Color'));
    
    for t = 1:nTimes % each time point
        subplot(1,nTimes,t);

        % Plot the Data 
        % -------------
        topoplot( meanData(:,t,c), chanlocs, 'style','map',... 'both'
            'maplimits', maplimits,... 'absmax',...
            'electrodes',p.Results.plotlabels, ... 'labels' / 'numpoint'
            'whitebk','on'...
            );
        
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
        if nTimes==1
            title([study(c).Condition ':  ' num2str(timePoints(:,t)') 'ms'], 'Interpreter', 'none');
        else
            if t==1
                if strcmpi(p.Results.t,'off')
                    suptitle([study(c).Condition]);
                else
                    suptitle(['DIFF:' conditions{1} '-' conditions{1}]);
                end
            end
            
            title([num2str(timePoints(:,t)') ' ms'], 'Interpreter', 'none');
        end
        
    end
    colormap(hot)
end

%% Save to R
if p.Results.R
    % Save the data in long format
    % ----------------------------
    nElec   = length(chanlocs);     % number of electrodes
    nTimes  = size(timePoints,2);   % number of time points
    nConds  = length(conditions);   % number of conditions
    
    % Get chanlocs Theta + Radius
    chanlocs    = struct2table(chanlocs);
    electrodes  = reshape(repmat(chanlocs{:,1},[1 nTimes nConds]),1,[])';
    theta       = reshape(repmat(chanlocs{:,8},[1 nTimes nConds]),1,[])';
    radius      = reshape(repmat(chanlocs{:,9},[1 nTimes nConds]),1,[])';
    
    % Repeat time points and labels
    conds       = reshape(repmat(conditions,[nElec*nTimes 1]),1,[])';
    if size(timePoints,1)==2
        for t = 1:nTimes
            temp{t} = num2str(timePoints(:,t)');
        end
        timePoints = temp;
        clear temp
    end
    times       = reshape(repmat(timePoints,[nElec nConds 1]),1,[])';
    amplitudes	= reshape(meanData,1,[])';
    
    % save to table:
    T   = table(conds,times,electrodes, theta, radius,amplitudes,....
        'VariableNames', {'Condition','TimePnt','Channel','Theta','Radius','Amp'});
    fn  = ['topoplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv
    
    
    % Make R code file
    % ----------------
    % save data + save r script
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