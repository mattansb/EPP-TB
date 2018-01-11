% This function find peak latencies and amplitudes, and is called internaly
% by measument functions.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
% See also epp_getamplitude & epp_getamplitude

%{
Change log:
-----------
25-11-2016  New function (written in MATLAB R2015a)
%}

function [amplitudes, latencies] = suppPeak(cutData, cutTime, direction, local)

for s = 1:size(cutData, 3) % for each subject
    for e = 1:size(cutData, 1) % for each electrode
        if local~=0 % if selected to find local peak
            ShiftRight = circshift(cutData(e,:,s), [0,1]); % shift data to the right
            ShiftLeft = circshift(cutData(e,:,s), [0,-1]); % shift data to the left
            
            if direction > 0
                cLocal = cutData(e,:,s) > ShiftRight & cutData(e,:,s) > ShiftLeft; % find points which are GREATER than both adjecent points
            elseif direction < 0
                cLocal = cutData(e,:,s) < ShiftRight & cutData(e,:,s) < ShiftLeft; % find points which are SMALLER than both adjecent points
            end
            
            if local > 1
                meanRight = tsmovavg(cutData(e,:,s), 's', local, 2);
                meanLeft = flip(tsmovavg(flip(cutData(e,:,s),2), 's', local, 2), 2);
                if direction > 0
                    cMean = cutData(e,:,s) > meanRight & cutData(e,:,s) > meanLeft; % find points which are GREATER than both adjecent means
                elseif direction < 0
                    cMean = cutData(e,:,s) < meanRight & cutData(e,:,s) < meanLeft; % find points which are SMALLER than both adjecent means
                end
                
                cLocal = cLocal & cMean;
            end
            cutData(e,[1:local end-local+1:end],s) = nan;
            
            cutData(e,~cLocal,s) = nan;
        end % end local
        
        if direction > 0
            [amplitudes(s,e), latencies(s,e)] = nanmax(cutData(e,:,s),[],2);
        elseif direction < 0
            [amplitudes(s,e), latencies(s,e)] = nanmin(cutData(e,:,s),[],2);
        end
        
        latencies(s,e) = cutTime(latencies(s,e));
    end
    latencies(isnan(amplitudes)) = nan;
end

end