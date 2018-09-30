% PURPOSE:  plot power and coherence (averaged across subjects and selected
%           electrodes), with or without error terms.
%
% FORMAT
% ------
% epp_plotTF(study,conditions,electrodes,varargin)
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
%           'scale'         - 'linear','log'. Y axis scale type (default:
%                             'linear'). 
%           'erspmaplimits' - [min max] values for the ersp color scale.
%           'itcmaplimits'  - [min max] values for the itc color scale.
%           'R'             - if true, plot data is saved into a csv file
%                             in the current directory + an R file with
%                             code to plot your data. (in R you can
%                             continue to format you plot - colors,
%                             annotations, etc.)
%
% See also epp_plotbutterfly, epp_plotgrands, epp_plottopo, epp_plottopoTF, epp_plotchannels
%
%
% Author: Mattan S. Ben Shachar, Rachel Rac & Michael Shmueli, BGU, Israel

%{
Change log:
-----------
14-05-2018  Improvment to exporting plot data
01-05-2018  Fixed bug that caused time-line to be flipped!
            Changed ersp colormap
15-03-2018  Fix and improvment to color limits
05-03-2018  Fix title printing
28-02-2018  ITC (abs) is now computed in function, to allow for the
            combination of conditions.
17-08-2017  New function (written in MATLAB R2015a)
%}

function epp_plotTF(study,conditions,electrodes,varargin)


%% Validate input data

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    addParameter(p,'scale', 'linear', @(x) any(strcmpi(x,{'linear','log'})))
    addParameter(p,'R', false, @islogical)
    addParameter(p,'erspmaplimits',nan,@isnumeric)
    addParameter(p,'itcmaplimits',nan,@isnumeric)
parse(p, study, conditions, electrodes, varargin{:}); % validate

%% Get only relevant conditions (in order!)

cInd    = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);
study   = study(cInd);

%% Get all Data
times   = study(1).timeLine;
frex    = study(1).freqs;
ytick   = round(logspace(log10(frex(1)),log10(frex(end)),10)*100)/100;

for c = 1:length(conditions)
    % Get ERSP
    ersp(c).data    = squeeze(mean(mean(study(c).ersp(electrodes,:,:,:),4),1))';
    ersp_max(c)     = max(max(ersp(c).data)); % find max point
    ersp_min(c)     = min(min(ersp(c).data)); % find min point
    
    % Get ITC
    if ~isreal(study(c).itc), study(c).itc = abs(study(c).itc); end
    itc(c).data     = squeeze(mean(mean(study(c).itc(electrodes,:,:,:),4),1))';
    itc_max(c)      = max(max(itc(c).data)); % find max point
    itc_min(c)      = min(min(itc(c).data)); % find min point
end

% Set ERSP color lims
if ~isnan(p.Results.erspmaplimits)
    ersp_range = p.Results.erspmaplimits;
else
    ersp_range = max([abs(min(ersp_min)) max(ersp_max)]);
    ersp_range = [-ersp_range ersp_range];
end

% Set ITC color lims
if ~isnan(p.Results.itcmaplimits)
    itc_range = p.Results.itcmaplimits;
else
    itc_range = [min(itc_min) max(itc_max)];
end
    
%% Plot

figure() % open new fig
for c = 1:length(conditions) % for each condition
    
    % Plot ERSP
    % ---------
    subplot(2,length(conditions),c);                                    % new subplot
    contourf(times,frex,ersp(c).data',40,'linecolor','none');
    set(gca,'ytick',ytick,'yscale',p.Results.scale)
    colormap(gca,f_makeColormap('wcbkryw'))
    caxis(ersp_range)
    
    if c == 1 % if this is the first plot
        xlabel('Time (ms)')
        ylabel('Frequency (Hz)')
        cl              = colorbar('west','AxisLocation','out');
        cl.Position(1)  = cl.Position(1)*0.4;
        cl.Position(3)  = cl.Position(3)*0.3;
        cl.Label.String = 'ERSP';

    end
    hold on
    title(conditions{c}, 'Interpreter', 'none')                     % add title
    plot([0 0],[frex(1) frex(end)],'Color', 'k', 'LineStyle', '-','LineWidth',1);
    plot([0 0],[frex(1) frex(end)],'Color', 'w', 'LineStyle', '-');
    
    
    % Plot ITC
    % --------
    subplot(2,length(conditions),length(conditions)+c)              % new subplot
    contourf(times,frex,itc(c).data',40,'linecolor','none')
    set(gca,'ytick',ytick,'yscale',p.Results.scale)
    colormap(gca,f_makeColormap('kryw'))
    caxis(itc_range)
    
    if c == 1 % if this is the first plot
        cl              = colorbar('west','AxisLocation','out');
        cl.Position(1)  = cl.Position(1)*0.4;
        cl.Position(3)  = cl.Position(3)*0.3;
        cl.Label.String = 'ITC';
    end
    hold on
    plot([0 0],[frex(1) frex(end)],'Color', 'k', 'LineStyle', '-','LineWidth',1);
    plot([0 0],[frex(1) frex(end)],'Color', 'w', 'LineStyle', '-');
end

%% Export to R?

if p.Results.R
    % Save the data in long format
    % ----------------------------
    save_data.Condition = {};
    save_data.Time      = [];
    save_data.Frequency = [];
    save_data.ersp      = [];
    save_data.itc       = [];

    nTmPnts = length(study(1).timeLine);    % length of time axis
    nFrex   = length(study(1).freqs);       % number of frequencies

    for c = 1:length(study) % for each condition
        temp_cond   = repmat({study(c).Condition},  [nTmPnts nFrex]);
        temp_time   = repmat(study(c).timeLine',    [1 nFrex]);
        temp_freq   = repmat(study(c).freqs,        [nTmPnts 1]);
        temp_ersp   = ersp(c).data;
        temp_itc    = itc(c).data;

        save_data.Condition = [save_data.Condition; temp_cond(:)];
        save_data.Time      = [save_data.Time;      temp_time(:)];
        save_data.Frequency = [save_data.Frequency; temp_freq(:)];
        save_data.ersp      = [save_data.ersp;      temp_ersp(:)];
        save_data.itc       = [save_data.itc;       temp_itc(:)];
    end

    T = struct2table(save_data);
    
    % Save to CSV
    % -----------
    fn  = ['waveletplot_' datestr(datetime, 'yyyymmdd_HHMMSS')];
    writetable(T,[fn '_data.csv'],'Delimiter',',','QuoteStrings',true) % save as csv
    
    % Make and save R code file
    % -------------------------
    Rpath = strrep(which(mfilename),[mfilename '.m'],'');
    fid  = fopen([Rpath '\epp_plotTF.R'],'rt');
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