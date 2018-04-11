% PURPOSE:  Download and install latest stable version of EPP-TB
%
% FORMAT
% ------
% download_install_EPPTB()
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
11-04-2018  New function (written in MATLAB R2017a)
%}
function download_install_EPPTB()

% Get instalation dir
save_path       = uigetdir([],'Select folder where EPP-TB will be installed');
save_zip        = [save_path '/EPP_TB.zip'];

if save_path==0
    error('Aborted download...')
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
    '\n\tWelcome to EPP-TB!'...
    '\n\tGet started by reading the <a href="https://github.com/mattansb/EPP-TB">README</a> file.'...
    '\n\t======================================='...
    '\n'...
    ])

end