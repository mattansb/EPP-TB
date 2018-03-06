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
% Optional arguments
% ------------------
% 'erp'         - if true, will compute mean-erps for each file.
% 'wavelet'     - if true, will compute ersp and itc for each file.
% 'savePath'    - char-vector of a folder in cd into which converted files
%                 will be saved. If the folder contains any '.eppf' files,
%                 these will be also be combined. (If the run colapses - as
%                 might happen with wavelet analysis - you can pick off
%                 from where it stopped.)
%                 If not specified, a folder with a random name will be
%                 created, and the name will be retuned.
% 'waveletVars' - a structure with parameters to be passed to
%                 suppWaveletConv3. If 'wavelet' is true, and 'waveletVars'
%                 is empty or not specified, a pop-up will asks for the
%                 parameters (which will also then be returned).
% 
%
%
% Author: Mattan S. Ben Shachar & Rachel Rac, BGU, Israel

%{
Change log:
-----------
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
parse(p, EEG_list, varargin{:}); % validate

% Wavelet parameters
% ------------------
if p.Results.wavelet
    if isempty(p.Results.waveletVars)
       [~, ~, ~, res] = inputgui(...
        'geometry', {[1.25 1 0.5] [1 1 0.75] [1 1 0.75] [1 0.5 1.25] [2.75] [2.75] [1 1 0.75] [1 1 0.75] [1 0.5 1.25]},...
        'geomvert', [1 1 1 1 1 1 1 1 1],...
        'title', 'wavelet',...
        'uilist', { ...
            {'Style', 'text', 'string', 'Frequencies and Cycles', 'fontweight', 'bold'  } {}...
            {'Style', 'pushbutton', 'string', 'help', 'callback', 'help suppWaveletConv3'} ...
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
end

% load eeglab
% -----------
eeglab redraw
close


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

for f = 1:nFiles
    %% Load
    temp_EEG = pop_loadset('filename',EEG_list{f});
    
    % condition name
    cond_parts = {temp_EEG.condition,temp_EEG.group};
    cond_parts = cond_parts(~cellfun(@isempty, cond_parts));
    
    output.Condition    = strjoin(cond_parts,'_');
    output.IDs          = table({temp_EEG.subject},temp_EEG.trials,'VariableNames',{'ID' 'nTrials'});
    
    %% ERP
    if p.Results.erp
        output.Data     = mean(temp_EEG.data,3);
        output.timeLine = temp_EEG.times;
    end
    
    %% Wavelet
    if p.Results.wavelet
        % to do!
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
    fname_parts = {output.Condition,output.IDs.ID{:}};
    fname_parts = fname_parts(~cellfun(@isempty, fname_parts));
    fname = [strjoin(fname_parts,'_') '.eppf'];
    
    if exist(fullfile(savePath, fname))==2
        output1 = output;
        load(fullfile(savePath, fname), '-mat')
        output(end+1) = output1;
    end
    
    save(fullfile(savePath, fname), 'output','-mat');
    
    clear output output1 fname fname_parts cond_parts temp_EEG
    
    
end

%% Combine all files

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
            study(c_ind).IDs             = [study(c_ind).IDs;output(c).IDs];
        end
    end
    clear output c_ind
end

end

function x = eval_res(x)
if ischar(x)
    eval(['x = [' x '];'])
else
    x = x;
end
end