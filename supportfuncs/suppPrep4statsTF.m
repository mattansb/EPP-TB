% This function prepares ERP data to be measured by the measurement
% functions - and is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
03-05-2018  Improvment to frequancy band selection
07-04-2018  New function (written in MATLAB R2015a)
%}
function [studyOut, TW, freqs_name] = suppPrep4statsTF(studyIn, conditions, electrodes, timeWindow, freqs, ave, jackknife)

%% Get only relevant conditions
cInd    = cellfun(@(x) find(ismember({studyIn(:).Condition}, x)), conditions);
studyIn = studyIn(cInd);

%% Find Closest time window
TW      = dsearchn(studyIn(1).timeLine',timeWindow(1)); % closest point to start of defined time window
TW(2)   = dsearchn(studyIn(1).timeLine',timeWindow(2)); % closest point to end of defined time window
TW      = TW(1):TW(2);

% Find Closest Frequencies
freqs_ind   = cell(1,size(freqs,1));
freqs_name  = cell(1,size(freqs,1));
for fr = 1:size(freqs,1)    
    F = find(studyIn(1).freqs >= freqs(fr,1) & studyIn(1).freqs <= freqs(fr,2));
    F = F([1 end]);
    
    freqs_ind{fr}   = F(1):F(2);
    freqs_name(fr)  = {[num2str(freqs(fr,1)) 'to' num2str(freqs(fr,2)) 'Hz']};
end

for c = 1:length(studyIn)
    %% Get freq-bands
    temp_ersp   = zeros(size(studyIn(c).ersp,1),size(freqs,1),size(studyIn(c).ersp,3),size(studyIn(c).ersp,4));
    temp_itc    = zeros(size(studyIn(c).itc,1),size(freqs,1),size(studyIn(c).itc,3),size(studyIn(c).itc,4));
    for fr = 1:length(freqs_ind)
        temp_ersp(:,fr,:,:) = mean(studyIn(c).ersp(:,freqs_ind{fr},:,:),2);
        temp_itc(:,fr,:,:)  = mean(studyIn(c).itc(:,freqs_ind{fr},:,:),2);
    end
    
    studyIn(c).ersp = temp_ersp;
    studyIn(c).itc  = temp_itc;
    
    %% Get only relevant electrodes
    studyIn(c).ersp = studyIn(c).ersp(electrodes,:,:,:);
    studyIn(c).itc  = studyIn(c).itc(electrodes,:,:,:);
    
    if ave
        studyIn(c).ersp = mean(studyIn(c).ersp,1);
        studyIn(c).itc  = mean(studyIn(c).itc,1);
    end
    
    %% Jackknife
    if jackknife
        studyIn(c).ersp = suppJackknife('in',studyIn(c).ersp,4);
        studyIn(c).itc  = suppJackknife('in',studyIn(c).itc,4);
    end
    
    %% Reduce dims   
    studyIn(c).origSize = size(studyIn(c).ersp);
    
    studyIn(c).ersp = permute(studyIn(c).ersp,[3 4 1 2]);
    studyIn(c).itc  = permute(studyIn(c).itc,[3 4 1 2]);
    
    studyIn(c).ersp = reshape(studyIn(c).ersp,size(studyIn(c).ersp,1),[]);
    studyIn(c).itc  = reshape(studyIn(c).itc,size(studyIn(c).itc,1),[]);
    
end

studyOut = studyIn;

end