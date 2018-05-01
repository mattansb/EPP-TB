% This function saves results obtained by the measurement functions - and
% is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
18-04-2018  Fix to bad chars in condition names (remove them)
19-01-2017  Support for saving data
19-12-2016  Support for Jackknife adjustment
25-11-2016  New function (written in MATLAB R2015a)
%}

function [results, study] = suppPrep4exportERP(measure,study,conditions,electrodes,pResults)

fprintf('\n\nSaving results..')



for c = 1:length(study) % for each condition
    %% reshape data
    study(c).measure = reshape(study(c).measure,study(c).origSize(3),[]);
    
    study(c).Data = reshape(study(c).Data,study(c).origSize(2),study(c).origSize(3),study(c).origSize(1));
    study(c).Data = permute(study(c).Data,[3 1 2]);
    
    %% Create results table for each condition
    
    % Correct for jackknifing:
    if pResults.jackknife
        study(c).measure = suppJackknife('out',study(c).measure,1);
    end
    
    % Make variable name:
    if pResults.average
        VariableNames = {[study(c).Condition '_ave']};
    else
        for e = 1:length(electrodes)
            VariableNames{e} = sprintf('%s_%d', conditions{c}, electrodes(e));
        end
    end
    % remove any bad chars from names
    VariableNames = regexp(VariableNames,'[0-9a-zA-Z_]+','match');
    VariableNames = cellfun(@(x) [x{:}],VariableNames,'UniformOutput',false);
    
    study(c).exp = array2table(study(c).measure, 'VariableNames', VariableNames);   % convert results to table
    study(c).exp = [study(c).IDs.ID,study(c).exp];                                  % merge results with IDS
end

%% Merge results tables to single table
results.(measure) = study(1).exp;
for c = 2:length(study)
    results.(measure) = outerjoin(results.(measure),study(c).exp,'MergeKeys',true);
end

%% Add measuemnt info
results.info = rmfield(pResults,'study');

%% Save to file?
if any(strcmpi(pResults.save, {'wide','long'}))
    fprintf('. Done!\n\n\n\nWriting to file..')
    
    save_data = results.(measure);
    
    if strcmpi(pResults.save, 'long')
        % get number of IDs and columns
        [nID, nCond] = size(save_data);
        nCond = nCond-1;
        
        % shape values
        save_vals = table2array(save_data(:,2:end));
        save_vals = reshape(save_vals,1,[])';
        
        % shape IDs
        save_ID = repmat(table2cell(save_data(:,1)),nCond,1);
        
        % shape condition column
        save_cond = save_data.Properties.VariableNames(2:end);
        save_cond = repmat(save_cond,nID,1);
        save_cond = reshape(save_cond,1,[])';
        
        % combine to table
        save_data = table(save_ID, save_cond, save_vals, 'VariableNames', {'ID','Condition',measure});
    end
        
    save_info = struct2table(results.info);
    
    fn  = ['erp_' measure '_' datestr(datetime, 'yyyymmdd_HHMMSS')]; % file name
    
    writetable(save_data,fn,'FileType','spreadsheet','Sheet',1) % write values
    writetable(save_info,fn,'FileType','spreadsheet','Sheet',2) % write info
end

fprintf('. Done!\n\n')

end