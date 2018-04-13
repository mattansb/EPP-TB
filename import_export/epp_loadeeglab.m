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
% EEG_list      - a cell-array of '.set' (full) file locations. If empty
%                 will skip, and will combine the files in 'savePath.
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
%                 suppWaveletConv3. If 'wavelet' is true, and 'waveletVars'
%                 is empty or not specified, a pop-up will asks for the
%                 parameters (which will also then be returned).
% 'combine'     - if true (defult) will combine all files in 'savePath'
%                 into study structure. Else, won't - will only conver
%                 eeglab files to '.eppf' files.
%
% See also epp_loaderplab, epp_loadegimat
%
%
% Author: Mattan S. Ben Shachar & Rachel Rac, BGU, Israel

%{
Change log:
-----------
13-04-2018  A million little bug fixes
11-04-2018  Added errors when subject ID
20-03-2018  Fix naming of subject files.
06-03-2018  New function (written in MATLAB R2017a)
%}
function [study, savePath, waveletVars] = epp_loadeeglab(EEG_list,varargin)

%% Validate and initiate
p = inputParser;
    addRequired(p,'EEG_list',@iscell);
    addParameter(p,'erp', false, @islogical)
    addParameter(p,'wavelet', false, @islogical)
    addParameter(p,'savePath', '', @ischar)
    addParameter(p,'waveletVars', [], @isstruct)
    addParameter(p,'combine', true, @islogical)
parse(p, EEG_list, varargin{:}); % validate

if ~(p.Results.wavelet || p.Results.erp || p.Results.combine)
    error('You seem to have called the function without asking it to make ERPs or Wavelets. Oops?')
end


% Wavelet parameters
% ------------------
if p.Results.wavelet
    if isempty(p.Results.waveletVars)
       [~, ~, ~, res] = inputgui(...
        'geometry', {[1.25 1 0.5] [1 1 0.75] [1 1 0.75] [1 0.5 1.25] [2.75] [2.75] [1 1 0.75] [1 1 0.75] [1 0.5 1.25]},...
        'geomvert', [1 1 1 1 1 1 1 1 1],...
        'title', 'Specify Wavelet Paramters',...
        'uilist', { ...
            {'Style', 'text', 'string', 'Frequencies and Cycles', 'fontweight', 'bold'  } {}...
            {'Style', 'pushbutton', 'string', 'help', 'callback', 'pophelp(''suppWaveletConv3'')'} ...
            {'Style', 'text', 'string', 'Frequency Range' }...
            {'Style', 'edit', 'string', '', 'tag' 'freqRange' }...
            {'Style', 'checkbox', 'string' 'log-space' 'value' 1 'tag' 'log' }...
            {'Style', 'text', 'string', 'Cycle Range' }...
            {'Style', 'edit', 'string', '', 'tag' 'cycleRange' } {} ...
            {'Style', 'text', 'string', 'N' }...
            {'Style', 'edit', 'string', '', 'tag' 'freqN' } {} ...
            {}...
            {'Style', 'text', 'string', 'Times', 'fontweight', 'bold'  }...
            {'Style', 'text', 'string', 'Baseline' }...
            {'Style', 'edit', 'string', '', 'tag' 'Baseline' } ...
            {'Style', 'checkbox', 'string' 'dB' 'value' 1 'tag' 'db' }...
            {'Style', 'text', 'string', 'Cut times' }...
            {'Style', 'edit', 'string', '', 'tag' 'Cut_times' } {}...
            {'Style', 'text', 'string', 'Downsample' }...
            {'Style', 'edit', 'string', '1', 'tag' 'Downsample' } {} ...
            } );


        waveletVars = structfun(@eval_res ,res,'UniformOutput',false);
    else
        waveletVars = p.Results.waveletVars;
    end
    
    res_wave = {waveletVars.freqRange,...
        waveletVars.freqN,...
        waveletVars.cycleRange,...
        waveletVars.Baseline,...
        waveletVars.Cut_times,...
        'log',waveletVars.log==1,...
        'dB',waveletVars.db==1,...
        'downsample',waveletVars.Downsample,...
        };
else
    waveletVars = p.Results.waveletVars;
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
        % Load and validate sub\group\condition
        % =====================================
        try
            temp_EEG = pop_loadset('filename',EEG_list{f});
        catch
            eeglab; close
            temp_EEG = pop_loadset('filename',EEG_list{f});
        end
        
        if isempty(temp_EEG.condition) && isempty(temp_EEG.group)
            error('EEG.condition and EEG.group are both missing (Need at least one of them).')
        elseif isempty(temp_EEG.subject)
            error('EEG.subject not speficied.')
        end
        
        cond_parts          = {temp_EEG.condition,temp_EEG.group};
        cond_parts          = cond_parts(~cellfun(@isempty, cond_parts));
        output.Condition    = strjoin(cond_parts,'_');
        output.IDs          = table({temp_EEG.subject},temp_EEG.trials,'VariableNames',{'ID' 'nTrials'});
        
        
        
        % ERP
        %====
        if p.Results.erp
            output.Data     = mean(temp_EEG.data,3);
            output.timeLine = temp_EEG.times;
        end

        % Wavelet
        % =======
        if p.Results.wavelet
            [power,itpc,frex,times] = suppWaveletConv3(temp_EEG,...
                res_wave{:},...
                'sound','off');

            output.ersp     = power;
            output.itc      = itpc;
            output.freqs    = frex;
            output.timeLine = times;

            clear power itpc frex times
        end

        %% Save
        fname = [output.IDs.ID{:} '.eppf'];

        if exist(fullfile(savePath, fname))==2
            output1 = output;
            load(fullfile(savePath, fname), '-mat')
            output(end+1) = output1;
        end

        save(fullfile(savePath, fname), 'output','-mat');

        clear output output1 fname fname_parts cond_parts temp_EEG
    end
end

%% Combine all files

if p.Results.combine
    % list files in savePath
    all_files   = dir([savePath '\\*.eppf']);
    nFiles      = length(all_files);

    study.Condition = '';
    for f = 1:nFiles
        % Load
        load(fullfile(savePath, all_files(f).name), '-mat')

        % Orgenize
        for c = 1:length(output)
            c_ind = find(strcmpi(output(c).Condition,{study.Condition}));
            if isempty(c_ind)
                if strcmpi(study(end).Condition,'')
                    study = output(c);
                else
                    study(end+1) = output(c);
                end
            else
                try study(c_ind).ersp(:,:,:,end+1)  = output(c).ersp;   end
                try study(c_ind).itc(:,:,:,end+1)   = output(c).itc;    end
                try study(c_ind).Data(:,:,end+1)    = output(c).Data;   end
                study(c_ind).IDs = [study(c_ind).IDs;output(c).IDs];
            end
        end
        clear output c_ind
    end
else
    % return an empty struct
    study = struct;
end


end

function x = eval_res(x)
if ischar(x)
    eval(['x = [' x '];'])
else
    x;
end
end