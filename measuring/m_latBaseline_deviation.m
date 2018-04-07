function res = m_latBaseline_deviation(data,timeWindow_ind,direction,baselineOS,times,SDcriterion)

try
    if direction < 0 % if max amp
        data = -1*data;
    end
    
    if baselineOS < 0
        baseline_ind = times >= baselineOS & times < 0;
    elseif baselineOS > 0
        baseline_ind = times <= baselineOS & times > 0;
    end

    Cr          = std(data(baseline_ind))*SDcriterion; % compute criterion = n*SD     
    cut_data    = data(timeWindow_ind);
    cut_times   = times(timeWindow_ind);
    
    ind = find(cut_data > Cr,1);   % find first time amp crosses criterion
    
    if ~isempty(ind)
        res = cut_times(ind);
    else
        res = nan;
    end
    
catch ME
    warning(ME.message)
    res = nan;
end

end