% This function saves results obtained by the measurement functions - and
% is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
% See also epp_getamplitude& & epp_getamplitude

%{
Change log:
-----------
07-04-2018  New function (written in MATLAB R2017a)
%}

function results = suppPrep4exportTF(measure,study,conditions,electrodes,freqs_name,pResults)

fprintf('\n\nSaving results..')

for c = 1:length(study) % for each condition
    study(c).measure = reshape(study(c).measure,...
        study(c).origSize(4),[]);
    
    %% Create results table for each condition
    
    % Correct for jackknifing:
    if pResults.jackknife
        warning('Not supported yet')
        % maybe add?
    end
    
    % Make variable name:
    VariableNames = {};
    for fr = 1:length(freqs_name)
        if pResults.average
            VariableNames{end+1} = {[study(c).Condition '_' freqs_name{fr} '_ave']};
        else
            for e = 1:length(electrodes)
                VariableNames{end+1} = sprintf('%s_%s_%d', conditions{c},freqs_name{fr}, electrodes(e));
            end
        end
    end
    
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