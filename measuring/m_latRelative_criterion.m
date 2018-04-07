function res = m_latRelative_criterion(data,timeWindow_ind,direction,local,times,percentage)

try
    [amp,lat] = m_Peak(data,timeWindow_ind,direction,local,times);
    
    if isnan(amp)
        res = nan;
        return
    end
    
    if direction < 0 % if max amp
        data    = -1*data;
        amp     = -1*amp;
    end
    
    Cr              = amp*percentage;       % compute criterion = %*peak_amp
    lat             = find(times==lat,1);   % convert laency back to index
    data(lat:end)   = nan;                  % remove later latencies
    
    ind = find(data < Cr,1,'last');   % find where amp is smaller than Criterion
    res = times(ind);
    
%     % Compute linear approximation:
%     cutTime = times(timeWindow_ind);
%     A   = data(ind+1)-data(ind);
%     a   = Cr - data(ind);
%     B   = cutTime(ind+1)-cutTime(ind);
%     res = cutTime(ind) + B*(a/A); % SAVE LATENCY!        
catch ME
    warning(ME.message)
    res = nan;
end

end