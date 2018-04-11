% PURPOSE:  Download and install latest stable version of EPP-TB
%
% FORMAT
% ------
% download_and_install_EPPTB()
%
% Process:
%   1.  Ask user to select a folder to install EPP-TB.
%   2.  Download the latest stable version of EPP-TB from GitHub.
%   3.  If a previous version is installed, will remove from Matlab path
%       (but wont delete the files.)
%   4.  Add the downloaded functions' paths to Matlab paths.
%   5.  Say hi (:
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
11-04-2018  New function (written in MATLAB R2017a)
%}
function download_and_install_EPPTB()

% Get instalation dir
save_path       = uigetdir([],'Select folder where EPP-TB will be installed');
save_zip        = [save_path '/EPP_TB.zip'];

if save_path==0
    error('Aborted download...')
end

% Rmoving old version
path_list   = strsplit(path,";")';
path_ind    = contains(path_list,'EPP-TB','IgnoreCase',true);
unin_txt    = '';
if any(path_ind)
    fprintf(['Uninstalling old version of EPP-TB\n'])
    EPP_path_list = path_list(path_ind);
    for p = 1:length(EPP_path_list)
        rmpath(EPP_path_list{p})
    end
    fprintf(['\b... Done.\n'])
    unin_txt = ['Previous version on EPP-TB was uninstalled, but files '...
        'may still be locally saved on your computer.)'];
end


% Download .zip file
fprintf(['Downloading .zip from GitHub\n'])
download_url    = 'https://api.github.com/repos/mattansb/EPP-TB/zipball';
outfilename     = websave(save_zip,download_url);
fprintf(['\b... Done.\n'])

% Unzip and delete zip
fprintf(['Unzipping folder\n'])
filenames = unzip(outfilename,save_path)';
eval(['delete ' outfilename])
fprintf(['\b... Done.\n'])

% Add to path
fprintf(['Adding to MATLAB path\n'])
[~,I] = min(cellfun(@length,filenames));
addpath(genpath(filenames{I}))
savepath
fprintf(['\b... Done.\n'])

% Say hi!
fprintf([...
    '\n\t======================================='...
    '\n\tWelcome to EPP-TB! (:'...
    '\n\tGet started by reading the <a href="https://github.com/mattansb/EPP-TB">README</a> file.'...
    '\n\t======================================='...
    '\n\n'...
    ])

if ~isempty(unin_txt)
    warning(unin_txt)
end

end