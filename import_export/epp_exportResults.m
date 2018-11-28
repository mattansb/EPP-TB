% PURPOSE:  measure ERP amplitudes.
%
% FORMAT
% ------
% epp_exportResults(results,save_format)
% 
%
% INPUTS
% ------
%           results     - a results structure returned by any epp_get*
%           save_format - 'long' / 'wide'; will save the results in the
%                         current directory in the specified format.
%
% See also epp_getAmp, epp_getLat, epp_getTF
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
28-11-2018  New function (written in MATLAB R2017b)
%}
function epp_exportResults(results,save_format)

fprintf('\nWriting to file..')

measure = fieldnames(results);
measure = measure {1};

save_data = results.(measure);

if strcmpi(save_format, 'long')
    save_data = stack(save_data,save_data.Properties.VariableNames(2:end),...
        'NewDataVariableName',measure,...
        'IndexVariableName','Condition');
    save_data.Condition = cellstr(save_data.Condition); % categorical to cellstr.
end

save_info = results.info;
try
    % for TF data
    save_info.freqs = cellfun(@(x) [num2str(x(1)) '-' num2str(x(2)) 'Hz'],...
    num2cell(save_info.freqs,2),...
    'UniformOutput',false)';
    fn  = ['wavelet_' measure '_' datestr(datetime, 'yyyymmdd_HHMMSS')]; % file name
catch
    fn  = ['erp_' measure '_' datestr(datetime, 'yyyymmdd_HHMMSS')]; % file name
end
save_info = cell2table(struct2cell(save_info)','VariableNames',fieldnames(save_info));

writetable(save_data,fn,'FileType','spreadsheet','Sheet',1) % write values
writetable(save_info,fn,'FileType','spreadsheet','Sheet',2) % write info

fprintf('. Done!\n\n')

end