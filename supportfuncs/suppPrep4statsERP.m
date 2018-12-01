% This function prepares ERP data to be measured by the measurement
% functions - and is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
07-04-2018  New function (written in MATLAB R2015a)
%}

function [studyOut, T] = suppPrep4statsERP(studyIn, conditions, electrodes, timeWindow, ave, jackknife, interpolate)

%% Get only relevant conditions
% cInd	= cellfun(@(x) find(ismember({studyIn(:).Condition}, x)), conditions);
cInd    = cellfun(@(x) find(strcmp(x,{studyIn(:).Condition})), conditions);
studyIn = studyIn(cInd);


%% Find Time Window
T    = dsearchn(studyIn(1).timeLine',timeWindow(1)); % closest point to start of defined time window
if length(timeWindow)==2
    T(2)    = dsearchn(studyIn(1).timeLine',timeWindow(2)); % closest point to end of defined time window
    T       = T(1):T(2);
end


for c = 1:length(studyIn)
    %% Get only relevant electrodes
    
    studyIn(c).Data = studyIn(c).Data(electrodes,:,:);
    if ave
        studyIn(c).Data = mean(studyIn(c).Data,1);        
    end
    
    %% Jackknife    
    if ~strcmp(jackknife,'off')
        if ~contains(jackknife,'unweighted')
            W = studyIn(c).IDs.nTrials;
        else
            W = 1;
        end
        
        [studyIn(c).Data, mean_dat] = f_jackknife('in',studyIn(c).Data,3,W);
        
        if ~contains(jackknife,'uncentered')
            studyIn(c).Data(:,:,end+1) = mean_dat;
        end
    end
    
    %% Interpolate
    if interpolate~=1 % if any need for interpolation
        samplingRate_2  = interpolate*1000/(studyIn(c).timeLine(2)-studyIn(c).timeLine(1));
        
        timeLine_1      = studyIn(c).timeLine;                                              % old time line
        timeLine_2      = timeLine_1(1):(1000/samplingRate_2):timeLine_1(end);              % new time line
        Data_1          = permute(studyIn(c).Data,[2,1,3]);                                 % old data
        Data_2          = permute(interp1(timeLine_1,Data_1,timeLine_2,'spline'),[2,1,3]);  % new data (interpolated!)
        
        studyIn(c).timeLine     = timeLine_2;       % save new time line
        studyIn(c).Data         = Data_2;           % save new data (interpolated!)
    end

    %% Reduce dims
    studyIn(c).origSize = size(studyIn(c).Data);
    studyIn(c).Data = permute(studyIn(c).Data,[2 3 1]);
    studyIn(c).Data = reshape(studyIn(c).Data,size(studyIn(c).Data,1),[]);
end

studyOut = studyIn;


end