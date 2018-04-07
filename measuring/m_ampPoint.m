function res = m_ampPoint(data,timeWindow_ind)

try
    res = data(timeWindow_ind);
catch ME
    warning(ME.message)
    res = nan;
end

end