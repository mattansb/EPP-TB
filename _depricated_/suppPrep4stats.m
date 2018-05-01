% This function prepares ERP data to be measured by the measurement
% functions - and is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
% See also epp_getamplitude& & epp_getamplitude

%{
Change log:
-----------
14-04-2018  Support new jackknife function
07-04-2018  Added samplingRate calculation
06-03-2018  Removed samplingRate & baseLine fields from study struct.
29-01-2017  Support for sampling interpolation
19-12-2016  Support for Jackknife adjustment
25-11-2016  New function (written in MATLAB R2015a)
%}

function studyOut = suppPrep4stats(studyIn, conditions, electrodes, timeWindow, ave, jackknife, interpolate)

%% Get only relevant conditions
% cInd	= cellfun(@(x) find(ismember({studyIn(:).Condition}, x)), conditions);
cInd    = cellfun(@(x) find(strcmp(x,{studyIn(:).Condition})), conditions);
studyIn = studyIn(cInd);

for c = 1:length(studyIn)
    %% Get only relevant electrodes
    if ave
        studyIn(c).Data = mean(studyIn(c).Data(electrodes,:,:),1);
    else
        studyIn(c).Data = studyIn(c).Data(electrodes,:,:);
    end
    
    %% Jackknife
    if jackknife
        studyIn(c).Data = suppJackknife('in',studyIn(c).Data,3);
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
    
    %% Find Time Window
    T(1)    = dsearchn(studyIn(c).timeLine',timeWindow(1)); % closest point to start of defined time window
    T(2)    = dsearchn(studyIn(c).timeLine',timeWindow(2)); % closest point to end of defined time window
    T(1)    = studyIn(c).timeLine(T(1));
    T(2)    = studyIn(c).timeLine(T(2));
    T       = studyIn(c).timeLine >= T(1) & studyIn(c).timeLine <= T(2);
        
    % Cut Data + Time
    studyIn(c).cutTime = studyIn(c).timeLine(T); % save the relevant time only
    studyIn(c).cutData = studyIn(c).Data(:,T,:); % save relevant data only
    
    clear T
end

studyOut = studyIn;


end