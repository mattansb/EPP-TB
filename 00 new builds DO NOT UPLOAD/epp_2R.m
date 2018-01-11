function T = epp_2R(study,chanlocs)

if isempty(chanlocs)
    chanlocs = 1:size(study(1).Data,1);
end

%% Shape data

for c = 1:length(study)
    % get sizes
    [nChan, nTime, nSub] = size(study(c).Data);
    
    % repeat electrods and times
    study(c).Condition = repmat({study(c).Condition},[nChan nTime nSub]);
    study(c).Channel = repmat(chanlocs',[1 nTime nSub]);
    study(c).Time = repmat(study(c).timeLine,[nChan 1 nSub]);
    study(c).Amplitude = study(c).Data;
    
    
    % get subject info
    for f = 1:size(study(c).IDs,2)
        vname = study(c).IDs.Properties.VariableNames{f};
        study(c).(vname) = permute(repmat(study(c).IDs{:,f},[1 nChan nTime]),[2 3 1]);        
    end
end
clear c f vname nChan nSub nTime chanlocs

study = rmfield(study,{'samplingRate','baseLine','IDs','timeLine','Data'});

for c = 1:length(study)
    fnames = fieldnames(study);
    for f = 1:length(fieldnames(study))
        study(c).(fnames{f}) = reshape(study(c).(fnames{f}),[],1);
    end
    
    % make table
    tempT = struct2table(study(c));
    tempT = tempT(:,[5 1:4 6:end]);
    
    % add to table
    if c == 1
        T = tempT;
    else
        T = [T;tempT];
    end
end
clear c f tempT fnames

%% Save
%{
aske where?
ask name?
%}


end