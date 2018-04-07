function res = m_latAbsolute_criterion(data,timeWindow_ind,direction,times,criterion)

try
    if direction < 0 % if max amp
        data        = -1*data;
        criterion   = abs(criterion);
    end
    
    cut_data    = data(timeWindow_ind);
    cut_times   = times(timeWindow_ind);
    
    ind = find(cut_data > criterion,1);   % find first time amp crosses criterion
    
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