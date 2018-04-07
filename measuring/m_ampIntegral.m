function res = m_ampIntegral(data,timeWindow_ind,samplingRate)

try
    res = sum(data(timeWindow_ind))/samplingRate; % compute integral mv/sec
catch ME
    warning(ME.message)
    res = nan;
end

end