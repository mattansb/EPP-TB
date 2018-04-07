function res = m_ampArea(data,timeWindow_ind,direction,samplingRate,boundary)

try
    data = data - boundary;
    if direction > 0      % only use positive values
        data(data < 0) = 0;
    elseif direction < 0  % only use negative values
        data(data > 0) = 0;
    end

    res = sum(abs(data(timeWindow_ind)))/samplingRate; % compute area mv/sec
catch ME
    warning(ME.message)
    res = nan;
end



end