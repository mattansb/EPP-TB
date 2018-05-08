% PURPOSE:  plot topos (averaged across subjects) using EEGLAB's topoplot
%           function.
%
% FORMAT
% ------
% epp_plottopoTF(study,chanlocs,conditions,timePoints,freqs,varargin)
%
%
% INPUTS
% ------
% Same as epp_plottopo, with the addition of:
% freqs         - matrix of frequencies, with each row containing a range
%                 of frequencies to group together (1st column is lower
%                 limit, 2nd column is upper limit of each range). e.g.
%                 freqs = [1 3; 4 15; 16 28];
%                 A given band is selected as so: [low <= freq < high]
%
% See also epp_plotbutterfly, epp_plotgrands, epp_plotTF, epp_plottopo
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
08-05-2018  Improvment to frequancy band selection
03-05-2018  Improvment to frequancy band selection
05-03-2018  New function (written in MATLAB R2017a)
%}

function epp_plottopoTF(study,chanlocs,conditions,timePoints,freqs,varargin)

%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'freqs',@(x) isnumeric(x) && size(x,2)==2);
parse(p, study, freqs); % validate
% all other parameters will be validated in epp_plottopo


%% Orgenize data

% Get only relevant conditions (in order!)
cInd  = cellfun(@(x) find(strcmp(x,{study(:).Condition})), conditions);
study = study(cInd);
clear cInd

% Find Closest Frequencies
nBands = size(freqs,1);

for fr = 1:nBands
    freqs_ind = find(study(1).freqs >= freqs(fr,1) & study(1).freqs < freqs(fr,2));
    freqs_ind = freqs_ind([1 end]);
    
    freqsRange(fr,:)    = study(1).freqs(freqs_ind);
    freqs_name(fr)      = {[num2str(freqs(fr,1)) 'to' num2str(freqs(fr,2)) 'Hz']};
end

%% Cut data
% for each frequency, average frequency and send to epp_plottopo as amp?
% add freqrange to conditon names, for better export to R from epp_plottopo

studyOut.Condition = '';
for fr = 1:nBands
    for c = 1:length(study)
        bool_freq = study(c).freqs >= freqsRange(fr,1) & study(c).freqs <= freqsRange(fr,2);
        
        studyOut(end+1).Condition = [study(c).Condition '_' freqs_name{fr}];
        studyOut(end).Data = squeeze(mean(study(c).ersp(:,bool_freq,:,:),2));
        studyOut(end).timeLine = study(c).timeLine;
    end
end

studyOut = studyOut(2:end);

%% Plot with epp_plottopo

epp_plottopo(studyOut,chanlocs,{studyOut.Condition},timePoints,varargin{:})


end