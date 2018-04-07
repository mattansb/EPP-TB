function res = m_tfITC(data,timeWindow_ind)

try
    res = mean(abs(data(timeWindow_ind))); % this is the itc data we want    
catch ME
    warning(ME.message)
    res = nan;
end

end