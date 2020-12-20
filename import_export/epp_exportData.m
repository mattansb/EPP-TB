% PURPOSE:  Export data to csv file
%
% Only support ERP data for now.
%
% FORMAT
% ------
% epp_exportData(study,varargin)
% 
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import.
%
%
% The available parameters are as follows:
%       'conditions'    - cell list of conditions to be saved. Must
%                         correspond to conditions in
%                         study(:).Condition.(e.g. {'freq', 'rare'}).
%                         Defaults to all conditions.
%       'electrodes'    - vector of electrodes to be saved after averaging
%                         (e.g. [87 85, 92]). Defaults to all electrodes.
%     	'average'       - average across electrodes before saving? (false
%                         by default).  
%       'timeWindow'    - A range of two time points, the data within will
%                         be saved. Defaults to all time points.
%
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
08-05-2020  New function (written in MATLAB R2017b)
%}
function epp_exportData(study,varargin)

%% Validate
if isfield(study,'ersp') || isfield(study,'itc')
    n_electrods = size(study(1).ersp, 1);
    f_range = study(1).freqs;
else
    n_electrods = size(study(1).Data, 1);
    f_range = nan;
end

p = inputParser;
    addRequired(p,'study',@isstruct);
    addParameter(p,'conditions', {study.Condition}, @iscellstr)
    addParameter(p,'electrodes', 1:n_electrods, @isnumeric)
    addParameter(p,'average', false, @islogical)
    addParameter(p,'timeWindow', study(1).timeLine([1, end]), @isnumeric)
    addParameter(p,'freqs', f_range, @isnumeric)
parse(p, study, varargin{:}); % validate


%% Clean data

% Conditions
% ----------
cInd    = cellfun(@(x) find(strcmp(x,{study(:).Condition})), p.Results.conditions);
study = study(cInd);

% electrodes
% ----------
if p.Results.average
    e_names = {'E_ave'};
else
    e_names = arrayfun(@(X) ['E' num2str(X)], p.Results.electrodes, 'UniformOutput', false);
end

for c = 1:length(study)
    if isfield(study,'ersp') || isfield(study,'itc')
        study(c).ersp   = study(c).ersp(p.Results.electrodes,:,:,:);
        study(c).itc    = abs(study(c).itc(p.Results.electrodes,:,:,:));
        
        if p.Results.average, study(c).ersp = mean(study(c).ersp, 1); end
        if p.Results.average, study(c).itc  = mean(study(c).itc, 1); end
    else
        study(c).Data = study(c).Data(p.Results.electrodes,:,:);
        if p.Results.average, study(c).Data = mean(study(c).Data, 1); end
    end
end


% timeWindow
% ----------
% closest point to start and end of defined time window
T = dsearchn(study(1).timeLine',p.Results.timeWindow'); 
T = T(1):T(2);
t_names = study(1).timeLine(T);

for c = 1:length(study)
    if isfield(study,'ersp') || isfield(study,'itc')
        study(c).ersp   = study(c).ersp(:,:,T,:);
        study(c).itc    = study(c).itc(:,:,T,:);
    else
        study(c).Data = study(c).Data(:,T,:);    
    end
end

% freqs
% ----------
% closest point to start and end of defined time window
if isfield(study,'ersp') || isfield(study,'itc')
    F = dsearchn(study(1).freqs',p.Results.freqs'); 
    f_names = study(1).freqs(F);
    
    study(c).ersp = study(c).ersp(:,F,:,:);
    study(c).itc = study(c).itc(:,F,:,:);
end

%% Tidy data

IDs = epp_extractIDs(study, {study.Condition});

outTable = table();

for c = 1:length(study)
    if isfield(study,'ersp') || isfield(study,'itc')
        error('Can''t do ersp/itc yet')
        % (condition * subject * time * frequency) * channels (ersp + itc)
        
        % id, condition, time and frequency
        
        % ersp and itc to array
        
        % append
    else
        % (condition * subject * time) * channels
    
        % id, condition, and time
        xid = repelem(study(c).IDs.ID, length(t_names), 1);
        xt  = repmat(t_names, 1, length(study(c).IDs.ID))';
        xc  = repmat({study(c).Condition}, 1, length(xt))';
        dat = table(xc,xid,xt,'VariableNames',{'Condition', 'ID', 'Time'});
    
        % Data to array
        eegdat = permute(study(c).Data, [1 2 3]);
        eegdat = reshape(eegdat,length(e_names),[]);
        eegdat = array2table(eegdat', 'VariableNames', e_names);
    
        % append
        dat = [dat, eegdat];
    
        outTable = [outTable; dat];
    end
end

outTable = outerjoin(outTable, IDs, 'Keys', {'ID','Condition'}, 'MergeKeys',true);

%% Save
% as csv

fn  = ['data_' datestr(datetime, 'yyyymmdd_HHMMSS') '.csv']; % file name
writetable(outTable,fn, 'QuoteStrings', true)

end