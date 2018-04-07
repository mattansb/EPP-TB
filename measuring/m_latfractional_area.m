function res = m_latfractional_area(data,timeWindow_ind,direction,boundary,times,percentage)

try
    data = data-boundary;

    % Rectify area
    if direction == 0 % all area
        data = abs(data); % make all amps positive!
    elseif direction > 0 % only positive area
        data(data < 0) = 0;
    elseif direction < 0 % only negative area
        data(data > 0)  = 0;
        data            = -1*data;
    end

    % turn amps in time window to cumulative area: 
    cut_data    = data(timeWindow_ind);
    cut_data    = cumsum(cut_data);
    cut_data    = cut_data/cut_data(end);
    cut_times   = times(timeWindow_ind);

    % compute criterion = %*area
    ind = find(cut_data > percentage,1);   % find first time amp crosses criterion
    
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