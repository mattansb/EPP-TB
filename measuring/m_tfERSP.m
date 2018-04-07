function res = m_tfERSP(data,timeWindow_ind)

try
    res = mean(data(timeWindow_ind)); % this is the ersp data we want
catch ME
    warning(ME.message)
    res = nan;
end

end