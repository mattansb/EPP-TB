% PURPOSE:  measure ERP amplitudes.
%
% FORMAT
% ------
% results = epp_getTF(measure,study,conditions, electrodes,timeWindow,freqs,varargin)
% 
%
% INPUTS
% ------
% measure       - 'ersp' / 'itc'
% study         - structure.
% conditions    - cell list of conditions to be plotted. Must correspond to
%                 conditions in study(:).Condition.(e.g. {'freq', 'rare'}).
% electrodes    - vector of electrodes to be plotted after averaging (e.g.
%                 [87 85, 92]).
% timeWindow    - two time points ([start end], in ms) within which the
%                 measument will be taken.
% freqs         - matrix of frequencies, with each row containing a range
%                 of frequencies to group together (1st column is lower
%                 limit, 2nd column is upper limit of each range). e.g.
%                 freqs = [1 3; 4 15; 16 28];
%
% The available parameters are as follows:
%           'average'       - average across electrodes before measuring?
%                             (false my default). 
%           'save'          - 'long' / 'wide'; will save the results in the
%                             current directory in the specified format.
%
%{
Change log:
-----------
28-02-2018  ITC (abs) is now computed in function, to allow for the
            combination of conditions.
17-08-2017  New function (written in MATLAB R2015a)

2DO
===
p.Results.plot??
%}
function results = epp_getTF(measure,study,conditions, electrodes,timeWindow,freqs,varargin)

%% Validate

p = inputParser;
    addRequired(p,'measure',@(x) any(strcmpi(x,{'ersp','itc'})));
    addRequired(p,'study',@isstruct);
    addRequired(p,'conditions',@iscellstr);
    addRequired(p,'electrodes',@(x) isvector(x) && isnumeric(x));
    addRequired(p,'timeWindow',@(x) isvector(x) && isnumeric(x) && length(x)==2);
    addRequired(p,'freqs',@(x) isnumeric(x) && size(x,2)==2);
    
    addParameter(p,'average',false,@islogical);
    addParameter(p,'plot',false,@islogical);
    addParameter(p,'save','no', @ischar);
parse(p, measure ,study, conditions, electrodes, timeWindow, freqs, varargin{:}); % validate


%% Orgenize Data Before Measuring

% Get only relevant conditions
cInd    = cellfun(@(x) find(ismember({study(:).Condition}, x)), conditions);
study   = study(cInd);

% Find Closest time window
timeWindow(1)   = dsearchn(study(1).timeLine',timeWindow(1)); % closest point to start of defined time window
timeWindow(2)   = dsearchn(study(1).timeLine',timeWindow(2)); % closest point to end of defined time window
timeWindow      = study(1).timeLine(timeWindow);
bool_time       = study(1).timeLine >= timeWindow(1) & study(1).timeLine <= timeWindow(2);

% Find Closest Frequencies
% fr = 1
for fr = 1:size(freqs,1)
    freqs1(fr,1) = dsearchn(study(1).freqs',freqs(fr,1)); % closest point to start of defined frequency window
    freqs1(fr,2) = dsearchn(study(1).freqs',freqs(fr,2)); % closest point to end of defined frequency window
    freqs1(fr,:)  = study(1).freqs(freqs1(fr,:));
    freqs_name(fr) = {[num2str(freqs(fr,1)) 'to' num2str(freqs(fr,2)) 'Hz']};
end

%% Get Measure

% c = 1
for c = 1:length(study)
    fprintf('\nCalculating for %s (condition %d of %d)..',study(c).Condition, c ,length(study))
%     fr = 1
    for fr = 1:size(freqs1,1)
        
        bool_freq = study(c).freqs >= freqs1(fr,1) & study(c).freqs <= freqs1(fr,2);
        
        switch measure
            case 'ersp'
                temp_data = study(c).ersp(electrodes,bool_freq,bool_time,:);    % this is the ersp data we want
            case 'itc'
                temp_data = study(c).itc(electrodes,bool_freq,bool_time,:);     % this is the itc data we want
                temp_data = abs(temp_data);
        end
        
        study(c).measure(:,fr,:) = squeeze(mean(mean(temp_data,2),3));    % compute mean: across electrodes, time window & freq window
    end
    
    
    
    if p.Results.average
        study(c).measure = permute(study(c).measure,[3 2 1]);         % reorder dimentions
        study(c).measure = mean(study(c).measure,3);
    else
        study(c).measure = permute(study(c).measure,[3 1 2]);         % reorder dimentions
    end
    fprintf('. Done!')
end

%% Save for export

fprintf('\n\nSaving results..')
for c = 1:length(study)
    
    for fr = 1:size(freqs,1)
        if p.Results.average
            measure_name{fr} = sprintf('%s_%s_ave', conditions{c}, freqs_name{fr});
        else
            for e = 1:length(electrodes)
                measure_name{fr,e} = sprintf('%s_%s_%de', conditions{c}, freqs_name{fr},electrodes(e));
            end
        end
    end
    
    if ~p.Results.average
        measure_name = reshape(measure_name',1,[]);
        study(c).measure = reshape(study(c).measure,size(study(c).measure,1),[]); 
    end
    
    study(c).exp = array2table(study(c).measure, 'VariableNames', measure_name);   % convert results to table
    study(c).exp = [study(c).IDs,study(c).exp];                                     % merge results with IDS
    
    clear measure_name
end

% Merge all conditions:
results.(measure) = study(1).exp;

for c = 2:length(study)
    results.(measure) = outerjoin(results.(measure),study(c).exp,'MergeKeys',true);
end

% Add measuemnt info
results.info = rmfield(p.Results,'study'); 

%% Save to file

if ~strcmpi(p.Results.save,'no')
    fprintf('. Done!\n\n\n\nWriting to file..')
    
    save_data = results.(measure);
    
    if strcmpi(p.Results.save,'long')
        % get number of IDs and columns
        [nID, nCond] = size(save_data);
        nCond = nCond-1;
        
        save_vals = table2array(save_data(:,2:end));
        save_vals = reshape(save_vals,1,[])';
        
        % shape IDs
        save_ID = repmat(table2cell(save_data(:,1)),nCond,1);
        
        % shape condition column
        save_cond = save_data.Properties.VariableNames(2:end);
        save_cond = repmat(save_cond,nID,1);
        save_cond = reshape(save_cond,1,[])';
        
        % combine to long table
        save_data = table(save_ID, save_cond, save_vals, 'VariableNames', {'ID','Condition',measure});
    end
    
    results.info.freqs = freqs_name;
    save_info = struct2table(results.info);
    
    
    fn  = ['wavelet_' measure '_' datestr(datetime, 'yyyymmdd_HHMMSS')]; % file name
    writetable(save_data,fn,'FileType','spreadsheet','Sheet',1) % write values
    writetable(save_info,fn,'FileType','spreadsheet','Sheet',2) % write info
end

fprintf('. Done!\n\n')

end