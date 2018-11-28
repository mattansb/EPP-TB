% This function saves results obtained by the measurement functions - and
% is called internaly by them when needed.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
07-07-2018  Minor change to output data type.
14-05-2018  Add 'ID' name to results table
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
    study(c).exp = [study(c).IDs(:,1),study(c).exp];                                % merge results with IDS
end

%% Merge results tables to single table
results.(measure) = study(1).exp;
for c = 2:length(study)
    results.(measure) = outerjoin(results.(measure),study(c).exp,'MergeKeys',true);
end

%% Add measuemnt info
results.info = rmfield(pResults,'study');
fprintf('. Done!\n\n')

%% Save to file?
if any(strcmpi(pResults.save, {'wide','long'}))
    epp_exportResults(results,pResults.save)
end

end