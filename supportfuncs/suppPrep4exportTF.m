% This function saves results obtained by the measurement functions - and
% is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
13-05-2018  Improvment to saving long data
02-05-2018  Improvemnt to saving measurement info that caused some bugs
18-04-2018  Fix to bad chars in condition names (remove them)
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
            VariableNames{end+1} = [study(c).Condition '_' freqs_name{fr} '_ave'];
        else
            for e = 1:length(electrodes)
                VariableNames{end+1} = sprintf('%s_%s_%d', conditions{c},freqs_name{fr}, electrodes(e));
            end
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
        save_data = stack(save_data,save_data.Properties.VariableNames(2:end),...
            'NewDataVariableName',measure,...
            'IndexVariableName','Condition');
    end
    
    save_info = results.info;
    save_info.freqs = cellfun(@(x) [num2str(x(1)) '-' num2str(x(2)) 'Hz'],...
        num2cell(save_info.freqs,2),...
        'UniformOutput',false)';
    save_info = cell2table(struct2cell(save_info)','VariableNames',fieldnames(save_info));
    
    fn  = ['wavelet_' measure '_' datestr(datetime, 'yyyymmdd_HHMMSS')]; % file name
    
    writetable(save_data,fn,'FileType','spreadsheet','Sheet',1) % write values
    writetable(save_info,fn,'FileType','spreadsheet','Sheet',2) % write info
end

fprintf('. Done!\n\n')

end