function [studyOut, TW, freqs_name] = suppPrep4statsTF(studyIn, conditions, electrodes, timeWindow, freqs, ave, jackknife)

%% Get only relevant conditions
cInd    = cellfun(@(x) find(ismember({studyIn(:).Condition}, x)), conditions);
studyIn = studyIn(cInd);

%% Find Closest time window
TW      = dsearchn(studyIn(1).timeLine',timeWindow(1)); % closest point to start of defined time window
TW(2)   = dsearchn(studyIn(1).timeLine',timeWindow(2)); % closest point to end of defined time window
TW      = TW(1):TW(2);

% Find Closest Frequencies
for fr = 1:size(freqs,1)    
    F(1) = dsearchn(studyIn(1).freqs',freqs(fr,1)); % closest point to start of defined frequency window
    F(2) = dsearchn(studyIn(1).freqs',freqs(fr,2)); % closest point to end of defined frequency window
    freqs_ind{fr}     = F(1):F(2);
    freqs_name(fr)      = {[num2str(freqs(fr,1)) 'to' num2str(freqs(fr,2)) 'Hz']};
end

for c = 1:length(studyIn)
    %% Get only relevant electrodes
    if ave
        studyIn(c).ersp = mean(studyIn(c).ersp(electrodes,:,:,:),1);
        studyIn(c).itc  = mean(studyIn(c).itc(electrodes,:,:,:),1);
    else
        studyIn(c).ersp = studyIn(c).ersp(electrodes,:,:,:);
        studyIn(c).itc  = studyIn(c).itc(electrodes,:,:,:);
    end
    
    %% Get freq-bands
    for fr = 1:length(freqs_ind)
        temp_ersp(:,fr,:,:) = mean(studyIn(c).ersp(:,freqs_ind{fr},:,:),2);
        temp_itc(:,fr,:,:)  = mean(studyIn(c).itc(:,freqs_ind{fr},:,:),2);
    end
    
    studyIn(c).ersp = temp_ersp;
    studyIn(c).itc  = temp_itc;
    
    %% Jackknife
    if jackknife
        warning('Not supported yet')
        % maybe add?
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