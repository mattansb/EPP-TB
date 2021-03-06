% PURPOSE:  import multiple eeglab files into an EPP-TB structure.
%
%
% FORMAT
% ------
% [study, savePath, waveletVars] = epp_loadeeglab(EEG_list,varargin)
%
%
% INPUTS
% ------
% EEG_list      One of the following:
%               1. A cell-array - (full) '.set' file locations. 
%               2. A string - name of a folder in which to select the
%                  '.set' files.
%               3. A logial [true] - will prompt a file selection window.
%               4. An empty cell - will skip to combine the files in 
%                  'savePath'.
%
% Condition names are taken from EEG.condition and EEG.group. Thus it is
% assumed that each file contains only a single condition.
% Subject IDs are taken from EEG.subject.
%
% 
% Optional arguments
% ------------------
% 'erp'         - if true, will compute mean-erps for each file.
% 'wavelet'     - if true, will compute ersp and itc for each file.
% 'savePath'    - char-vector of a folder in cd into which converted files
%                 will be saved. If the folder contains any '.eppf' files,
%                 these will be also be combined. (If the run crashes - as
%                 might happen with wavelet analysis - you can pick off
%                 from where it stopped.)
%                 If not specified, a folder with a random name will be
%                 created, and the name will be retuned.
% 'waveletVars' - a structure with parameters to be passed to
%                 f_WaveletConv. If 'wavelet' is true, and 'waveletVars'
%                 is empty or not specified, a pop-up will asks for the
%                 parameters (which will also then be returned).
% 'combine'     - if true (defult) will combine all files in 'savePath'
%                 into study structure. Else, won't - will only convert
%                 eeglab files to '.eppf' files.
% 'split'       - If true, will split each condition to '*_splitOdd' and
%                 '*_splitEven'. This is useful for measuring split half
%                 relability.
%
% See also epp_loaderplab, epp_loadegimat
%
%
% Author: Mattan S. Ben Shachar & Rachel Rac, BGU, Israel

%{
Change log:
-----------
21-05-2020  Added informative error when eeglab not attached.
05-05-2020  Added more options for file imports.
            Show warning when the same condition is being loaded more than
            once.
18-11-2018  Better error when eeglab not loaded
06-07-2018  Minor printing adjustment
10-05-2018  Support for new baseline correction methods in f_WaveletConv
25-04-2018  Fix when combining
16-04-2018  Improvments to speed when combining data
15-04-2018  Reduced printing from eeglab pop_*
13-04-2018  A million little bug fixes
11-04-2018  Added errors when subject ID
20-03-2018  Fix naming of subject files.
06-03-2018  New function (written in MATLAB R2017a)
%}
function [study, savePath, waveletVars] = epp_loadeeglab(EEG_list,varargin)

%% Validate and initiate
p = inputParser;
    addParameter(p,'erp', false, @islogical)
    addParameter(p,'wavelet', false, @islogical)
    addParameter(p,'savePath', '', @ischar)
    addParameter(p,'waveletVars', [], @isstruct)
    addParameter(p,'combine', true, @islogical)
    addParameter(p,'split', false, @islogical)
parse(p, varargin{:}); % validate

if ~(p.Results.wavelet || p.Results.erp || p.Results.combine)
    error('You seem to have called the function without asking it to make ERPs or Wavelets. Oops?')
end

if (p.Results.wavelet || p.Results.erp) && exist('pop_loadset','file')~=2
    error('eeglab must be loaded for this function to work...')
end

if ~iscell(EEG_list)
    if islogical(EEG_list) && EEG_list
        [file,path,indx] = uigetfile('*.set', 'Select EEGLAB set files','MultiSelect','on');
        
        if ~indx, error('No files selected.'); end
        
        EEG_list = fullfile(path,file);
    elseif ischar(EEG_list)
        listing = dir(EEG_list);
        listing = fullfile({listing.folder},{listing.name});
        
        ix = regexp(listing, '\.set$','start','forceCellOutput');
        ix = ~cellfun(@isempty, ix);
        
        EEG_list = listing(ix);
    else
        error('Wrong input for EEG_list. See help.')
    end
end


% Wavelet parameters
% ------------------
if p.Results.wavelet
    if isempty(p.Results.waveletVars)
       baseline_methods = {'dB','normalize','standardize'};
       [~, ~, ~, res] = inputgui(...
        'geometry', {[1.25 1 0.5] [1 1 0.75] [1 1 0.75] [1 0.5 1.25] [2.75] [2.75] [1 1 0.75] [1 1 0.75] [1 0.5 1.25]},...
        'geomvert', [1 1 1 1 1 1 1 1 1],...
        'title', 'Specify Wavelet Paramters',...
        'uilist', { ...
            {'Style', 'text', 'string', 'Frequencies and Cycles', 'fontweight', 'bold'  } {}...
            {'Style', 'pushbutton', 'string', 'help', 'callback', 'pophelp(''f_WaveletConv'')'} ...
            {'Style', 'text', 'string', 'Frequency Range' }...
            {'Style', 'edit', 'string', '', 'tag' 'freqRange' }...
            {'Style', 'checkbox', 'string' 'log-space' 'value' 1 'tag' 'log' }...
            {'Style', 'text', 'string', 'Cycles per Wavelet' }...
            {'Style', 'edit', 'string', '', 'tag' 'cycles' } ...
            {'Style', 'text', 'string', '(constant or range)' }...
            {'Style', 'text', 'string', 'N' }...
            {'Style', 'edit', 'string', '', 'tag' 'freqN' } {} ...
            {}...
            {'Style', 'text', 'string', 'Times', 'fontweight', 'bold'  }...
            {'Style', 'text', 'string', 'Baseline' }...
            {'Style', 'edit', 'string', '', 'tag' 'Baseline_time' } ...
            {'Style', 'popup', 'string' baseline_methods 'tag' 'baseline_method' }...
            {'Style', 'text', 'string', 'Cut times' }...
            {'Style', 'edit', 'string', '', 'tag' 'Cut_times' } {}...
            {'Style', 'text', 'string', 'Downsample' }...
            {'Style', 'edit', 'string', '1', 'tag' 'Downsample' } {} ...
            } );
        
        waveletVars = structfun(@eval_res ,res,'UniformOutput',false);
        waveletVars.baseline_method = baseline_methods{waveletVars.baseline_method};
    else
        waveletVars = p.Results.waveletVars;
    end
    
    res_wave = {waveletVars.freqRange,...
        waveletVars.freqN,...
        waveletVars.cycles,...
        waveletVars.Baseline_time,...
        waveletVars.Cut_times,...
        'log',waveletVars.log==1,...
        'baseline',waveletVars.baseline_method,...
        'downsample',waveletVars.Downsample,...
        };
end

%% Set up temp dir

% folder name
savePath = p.Results.savePath;
if isempty(savePath)
    symbols     = ['a':'z' 'A':'Z' '0':'9'];
    savePath    = symbols(randi(numel(symbols),[1 10]));
end

% create if doesnt exist
if ~exist(savePath,'dir')
    status = mkdir(savePath);
end

%% Run on all files

nFiles = length(EEG_list);

if p.Results.wavelet || p.Results.erp    
    for f = 1:nFiles
        fprintf('Processing file %d of %d\n',f,nFiles)
        % Load and validate sub\group\condition
        % =====================================
        try
            evalc('temp_EEG = pop_loadset(''filename'',EEG_list{f});');
        catch err
            error(err.message)
        end
        
        if isempty(temp_EEG.condition) && isempty(temp_EEG.group)
            error('EEG.condition and EEG.group are both missing (Need at least one of them).')
        elseif isempty(temp_EEG.subject)
            error('EEG.subject not speficied.')
        end
        
        cond_parts          = {temp_EEG.condition, temp_EEG.group};
        cond_parts          = cond_parts(~cellfun(@isempty, cond_parts));        
        
        if p.Results.split
            n_trial = temp_EEG.trials;
            i_odd   = 1:2:n_trial;
            i_even  = 2:2:n_trial;
            
            output(1).Condition    = [strjoin(cond_parts,'_') '_splitOdd'];
            output(2).Condition    = [strjoin(cond_parts,'_') '_splitEven'];
            
            output(1).IDs          = table({temp_EEG.subject},length(i_odd),'VariableNames',{'ID' 'nTrials'});
            output(2).IDs          = table({temp_EEG.subject},length(i_even),'VariableNames',{'ID' 'nTrials'});
            
            evalc('go_EEG = pop_select( temp_EEG, ''trial'',i_odd);');
            evalc('go_EEG(2) = pop_select( temp_EEG, ''trial'',i_even);');
        else
            output.Condition    = strjoin(cond_parts,'_');
            output.IDs          = table({temp_EEG.subject},temp_EEG.trials,'VariableNames',{'ID' 'nTrials'});
            
            go_EEG = temp_EEG;
        end
            
        for sp = 1:length(go_EEG)
            % ERP
            % ===
            if p.Results.erp
                output(sp).Data     = mean(go_EEG(sp).data,3);
                output(sp).timeLine = go_EEG(sp).times;
            end

            % Wavelet
            % =======
            if p.Results.wavelet
                [power,itpc,frex,times] = f_WaveletConv(go_EEG(sp),res_wave{:});
    
                output(sp).ersp     = power;
                output(sp).itc      = itpc;
                output(sp).freqs    = frex;
                output(sp).timeLine = times;
    
                clear power itpc frex times
            end
        end
            
        

        %% Save
        fname = [output(1).IDs.ID{:} '.eppf'];

        if exist(fullfile(savePath, fname))==2
            output1 = output;
            load(fullfile(savePath, fname), '-mat')
            match_cond = cellfun(@(cond) any(strcmp(cond,{output.Condition})), {output1.Condition});
            if any(match_cond)
                error(sprintf(['Condition ' output1(match_cond).Condition ' appears more than once for ID ' num2str(temp_EEG.subject)]))
            end
            output = [output, output1];
        end

        save(fullfile(savePath, fname), 'output','-mat');

        clear output output1 fname fname_parts cond_parts temp_EEG go_EEG
    end
end

%% Combine all files

study = struct;

if p.Results.combine
    % list files in savePath
    all_files   = dir([savePath '\\*.eppf']);
    nFiles      = length(all_files);

    study.Condition = '';
    % Load and append data
    % --------------------
    fprintf('\n\nLoading file ');
    nb = 0;
    for f = 1:nFiles
        fprintf(repmat('\b',1,nb));
        nb = fprintf('%d of %d\n',f,nFiles);
        % Load
        load(fullfile(savePath, all_files(f).name), '-mat')

        % Orgenize
        for c = 1:length(output)
            c_ind = find(strcmpi(output(c).Condition,{study.Condition}));
            if isempty(c_ind) % if there is no such condition
                if isempty(study(end).Condition) % if this is the first condition
                    study = output(c);
                else
                    study(end+1) = output(c);
                end
                
                % build data
                if isfield(study,'Data')
                    study(end).Data     = cell([1 nFiles]);
                    study(end).Data{1}  = output(c).Data;  
                    
                end
                
                % build ersp
                if isfield(study,'ersp')
                    study(end).ersp     = cell([1 nFiles]);
                    study(end).ersp{1}  = output(c).ersp;  
                end
                
                % build itc
                if isfield(study,'itc')
                    study(end).itc      = cell([1 nFiles]);
                    study(end).itc{1}   = output(c).itc;  
                end
                    
                
            else % condition exists
                study(c_ind).IDs = [study(c_ind).IDs;output(c).IDs];
                
                % append data
                if isfield(study,'Data')
                    next_ind = find(cellfun(@isempty, study(c_ind).Data),1);
                    study(c_ind).Data{next_ind} = output(c).Data;  
                end
                
                % append ersp
                if isfield(study,'ersp')
                    next_ind = find(cellfun(@isempty, study(c_ind).ersp),1);
                    study(c_ind).ersp{next_ind} = output(c).ersp;  
                end
                
                % append itc
                if isfield(study,'itc')
                    next_ind = find(cellfun(@isempty, study(c_ind).itc),1);
                    study(c_ind).itc{next_ind} = output(c).itc;  
                end                
            end
        end
        clear output c_ind
    end
    fprintf('\b... done.');
    
    % Collaps Data/ersp/itc
    % ---------------------
    fprintf('\nCleaning up (this can take a while)..')
    for c = 1:length(study)
        % Collaps data
        if isfield(study,'Data')
            next_ind = cellfun(@(X) ~isempty(X), study(c).Data);
            study(c).Data = cat(3,study(c).Data{next_ind});
        end

        % Collaps ersp
        if isfield(study,'ersp')
            next_ind = cellfun(@(X) ~isempty(X), study(c).ersp);
            study(c).ersp = cat(4,study(c).ersp{next_ind});
        end

        % Collaps itc
        if isfield(study,'itc')
            next_ind = cellfun(@(X) ~isempty(X), study(c).itc);
            study(c).itc = cat(4,study(c).itc{next_ind});
        end
    end
    fprintf('. Done!\n')
end


end

function x = eval_res(x)
if ischar(x), eval(['x = [' x '];']); end
end