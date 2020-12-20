% PURPOSE:  Convert TF data to "ERP" like data structure.
%
% FORMAT
% ------
% study = epp_reshapeTF(study, freqs, type)
%
%
% INPUTS
% ------
% study         - structure built by epp_load OR epp_erplab_import with
%                 ersp and itc data.
% freqs         - matrix of frequencies, with each row containing a range
%                 of frequencies to group together (1st column is lower
%                 limit, 2nd column is upper limit of each range). e.g.
%                 freqs = [1 3; 4 15; 16 28];  
%                 A given band is selected as so: [low <= freq <= high]
% type'         - 'erps' or 'itc'.
%
% Author: Mattan S. Ben Shachar, BGU, Israel

%{
Change log:
-----------
21-05-2020  New function (written in MATLAB R2017b)
%}
function studyOut = epp_reshapeTF(studyIn, freqs, type)
%% Validate

p = inputParser;
    addRequired(p,'study',@isstruct);
    addRequired(p,'freqs',@(x) isnumeric(x) && size(x,2)==2);
    addRequired(p,'type',@(x) strcmpi(x,'ersp') || strcmpi(x,'itc'))
parse(p, studyIn, freqs, type); % validate
% all other parameters will be validated in epp_plottopo


%% Orgenize data

% Find Closest Frequencies
nBands = size(freqs,1);
freqs_name = cell(nBands, 1);
freqsRange = nan(nBands, 2);

for fr = 1:nBands
    freqs_ind = find(studyIn(1).freqs >= freqs(fr,1) & studyIn(1).freqs <= freqs(fr,2));
    freqs_ind = freqs_ind([1 end]);
    
    freqsRange(fr,:)    = studyIn(1).freqs(freqs_ind);
    freqs_name(fr)      = {[num2str(freqs(fr,1)) 'to' num2str(freqs(fr,2)) 'Hz']};
end

%% Cut data
% for each frequency, average frequency and send to epp_plottopo as amp?
% add freqrange to conditon names, for better export to R from epp_plottopo

studyOut.Condition = '';
for fr = 1:nBands
    for c = 1:length(studyIn)
        bool_freq = studyIn(c).freqs >= freqsRange(fr,1) & studyIn(c).freqs <= freqsRange(fr,2);
        
        studyOut(end+1).Condition   = [studyIn(c).Condition '_' freqs_name{fr}];
        if strcmpi(type, 'ersp')
            studyOut(end).Data          = squeeze(mean(studyIn(c).ersp(:,bool_freq,:,:),2));
        else 
            if ~isreal(studyIn(c).itc), studyIn(c).itc = abs(studyIn(c).itc); end
            studyOut(end).Data          = squeeze(mean(studyIn(c).itc(:,bool_freq,:,:),2));
        end
        studyOut(end).timeLine      = studyIn(c).timeLine;
        studyOut(end).IDs           = studyIn(c).IDs;
    end
end

studyOut = studyOut(2:end);

end