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
02-05-2018  Improvemnt to saving measurement info
18-04-2018  Fix to bad chars in condition names (remove them)
19-01-2017  Support for saving data
19-12-2016  Support for Jackknife adjustment
25-11-2016  New function (written in MATLAB R2015a)
%}

function [results, study] = suppPrep4exportERP(measure,study,conditions,electrodes,pResults)

fprintf('Saving results..')



for c = 1:length(study) % for each condition
    %% reshape data
    study(c).measure = reshape(study(c).measure,study(c).origSize(3),[]);
    
    study(c).Data = reshape(study(c).Data,study(c).origSize(2),study(c).origSize(3),study(c).origSize(1));
    study(c).Data = permute(study(c).Data,[3 1 2]);
    
    %% Create results table for each condition
    
    % Correct for jackknifing:    
    if ~strcmp(pResults.jackknife,'off')
        
        
        
        if ~contains(pResults.jackknife,'unweighted')
            W = study(c).IDs.nTrials;    
        else
            W = 1;
        end
        
        if ~contains(pResults.jackknife,'uncentered')
            WM = study(c).measure(end);
            study(c).measure(end)   = [];
            study(c).Data(:,:,end)  = [];
        else
            WM = [];
        end
        
        study(c).measure = f_jackknife('out',study(c).measure,1,W,WM);
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
    study(c).exp = [study(c).IDs(:,1),study(c).exp];                                % merge results with IDS
end

%% Merge results tables to single table
results.(measure) = study(1).exp;
for c = 2:length(study)
    results.(measure) = outerjoin(results.(measure),study(c).exp,'MergeKeys',true);
end

%% Add measuemnt info
results.info = rmfield(pResults,'study');
fprintf('. Done!\n')

%% Save to file?
if any(strcmpi(pResults.save, {'wide','long'}))
    epp_exportResults(results,pResults.save)
end


end