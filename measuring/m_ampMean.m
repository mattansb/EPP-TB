function res = m_ampMean(data,timeWindow_ind)

try
    res = mean(data(timeWindow_ind)); % compute mean
catch ME
    warning(ME.message)
    res = nan;
end

end