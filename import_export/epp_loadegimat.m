% PURPOSE:  converts EGI mat file to a workable structure with the epp-tb
%           package.
%
% FORMAT
% ------
% study = epp_loadegimat(type,nchan,baseLine)
%
% 
%
% INPUTS
% ------
% type      - 'multi': RECOMNDED.
%                       will open a file selection window for selecting
%                       sevral indevidual .mat files (exported from NS).
%                       All files will be loaded into the STUDY structure,
%                       with subject IDs taken from the file names.
%             'combined': NOT RECOMNDED.
%                       will open a file selection window for selecting a
%                       SINGLE .mat file made from a combined file. NO ID
%                       LIST IS CONSTRUCTED, AND ONE MUST BE MADE MANUELY!
% nchan     - number of channels used in recording (e.g. 128).
% baseLine  - length of base line (negative or positive)
%
% OUTPUT
% ------
% study     - a structure orgenized by condition, with the following
%             feilds:
%                   1. Condition - name of condition (string)
%                   2. Data - ERP data, orgenized in 3D matrix
%                      (electrode*time*subject).
%                   3. samplingRate - sampling rate in Hz (e.g. 250)
%                   4. baseLine - same as input.
%                   5. timeLine - vector of time points corresponding to
%                      the time dimention in study.Data matrix.
%                   6. IDs - table containing IDs of subject in the
%                      condition. (ONLY IF USING THE 'MULTI' OPTION)
%
% See also epp_loaderplab, epp_loadeeglab
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
06-03-2018  Removed samplingRate & baseLine fields from study struct.
07-05-2017  Minor fix for importing data that was exported badly from NS
22-01-2017  Minor fix for selecting data
26-11-2016  Fix ID names, so saved without file extention
25-11-2016  New function (written in MATLAB R2015a)
%}

function study = epp_loadegimat(type,nchan,baseLine)

%% Validate
p = inputParser;
    addRequired(p,'type',@ischar);
    addRequired(p,'nchan',@isscalar);
    addRequired(p,'baseLine',@(x) isscalar(x) && isnumeric(x));
parse(p,type,nchan,baseLine);

switch lower(type)
%% For Combined File
case 'combined'
[File.Name,File.Path,~] = uigetfile('*.mat','File Selector','MultiSelect','off');

% Get list of condition names:
S           = whos('-file',fullfile(File.Path, File.Name));
CondVars    = cellfun(@(x) (x(1)==nchan || x(1)==nchan+1) && x(2)>1, {S.size})'; % only fields that have data with enough "electrodes"
CondVars    = {S(CondVars).name};

% Pre load the mat-file:
matF = matfile(fullfile(File.Path, File.Name));

for c = 1:length(CondVars) % for each condition
    ERPs(c).Condition       = CondVars{c};          % save condition name
    ERPs(c).Data            = matF.(CondVars{c});   % load data
            
    % Calculate time-line baed on base-line and sampling rate:
    Hz = 1000/(matF.samplingRate);
    if baseLine<0
        ERPs(c).timeLine    = baseLine:Hz:Hz*size(ERPs(c).Data,2)+baseLine-1; % this is better than using linspace(bl,x2,n)
    elseif baseLine>0
        ERPs(c).timeLine    = (-1)*(Hz*size(ERPs(c).Data,2)-baseLine):Hz:baseLine-1;
    end
    
    warning('When importing a combined file, subject IDs are not appended!')
end


%% For Individual Subject Files
case 'multi'
[File.Name,File.Path,~] = uigetfile('*.mat','File Selector','MultiSelect','on');

ERPs.Condition = '';

for s = 1:length(File.Name)
    fprintf('Loading file %d (of %d)...',s ,length(File.Name))
    
    % Get list of condition names:
    S           = whos('-file',fullfile(File.Path, File.Name{s}));
    CondVars    = cellfun(@(x) (x(1)==nchan || x(1)==nchan+1) && x(2)>1, {S.size})'; % only fields that have data with enough "electrodes"
    CondVars    = {S(CondVars).name};
    
    % Pre load the mat-file:
    matF = matfile(fullfile(File.Path, File.Name{s}));
    
    for c = 1:length(CondVars) % for each condition
        if size(matF.(CondVars{c}),3)>1
            temp_data = mean(matF.(CondVars{c}),3);
        else
            temp_data = matF.(CondVars{c});
        end
        
        if ~any(strcmpi(CondVars{c},{ERPs(:).Condition})) % if condition hasn't been previously loaded:
            ERPs(end+1).Condition   = CondVars{c};          % add new condition (name)
            ERPs(end).Data          = temp_data;   % load data
                       
            % Calculate time-line baed on base-line and sampling rate:
            Hz = 1000/(matF.samplingRate);
            if baseLine<0
                ERPs(end).timeLine    = baseLine:Hz:Hz*size(ERPs(end).Data,2)+baseLine-1; % this is better than using linspace(bl,x2,n)
            elseif baseLine>0
                ERPs(end).timeLine    = (-1)*(Hz*size(ERPs(end).Data,2)-baseLine):Hz:baseLine-1;
            end
            
            [~, f]                      = fileparts(File.Name{s});  % get file name without extention
            ERPs(end).IDs               = table({f},'VariableNames',{'ID'}); % save ID name
        else % if condition HAS been previously loaded:
            ind = strcmpi(CondVars{c},{ERPs(:).Condition});         % find where in the struct the condition was saved
            
            ERPs(ind).Data(:,:,end+1)	= temp_data;       % append new data to end of data
            
            [~, f]                      = fileparts(File.Name{s});  % get file name without extention
            ERPs(ind).IDs{end+1,1}      = {f};                      % append new IDto end of ID list
        end
    end
    fprintf('Done!\n')
end

ERPs    = ERPs(2:end);  % remove first (blank) row.

end

study   = ERPs;         % save for export

fprintf('Done!\n')

end