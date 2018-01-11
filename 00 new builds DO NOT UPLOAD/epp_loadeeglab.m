% must have eeglab loaded
function [study, strudy_rel] = epp_loadeeglab(inlist,type,varargin)

%{
varargin
save? pop out names
if type wavelet / both - pop up input to wavelet_conv2 - or... input as a
cell array?

%}

%%
% if isempty(inlist)
%     % open file selector?
% end

%% Build Struct

study.Condition = '';
for f = 1:length(inlist)
    % Load set file
    EEG = pop_loadset('filename',inlist{f});
    
    % initiate
    ID      = EEG.subject;
    cond    = EEG.condition;
    ntrials = EEG.trials;
    
    if strcmpi(type,{'ERP','both'})
        data = mean(EEG.data,3);
    end
    
    if strcmpi(type,{'WAVELET','both'})
        % NEED TO BUILD DATA
        [ERSP,ITC,frex,~] = wavelet_conv2(); % INPUT????????????????????
    end
    
    c_ind = strcmpi(cond,{study(:).Condition});
    
    if isempty(c_ind)
        study(end+1).Condition  = cond;
        study(end).samplingRate = EEG.srate;
        study(end).timeLine     = EEG.times;
        study(end).ID           = table(ID,ntrials,'VariableNames',{'ID' 'nTrials'});
        
        switch lower(type)
            case 'erp'
                study(end).Data     = ERP;
            case 'wavelet'
                study(end).Freqs    = frex;
                study(end).ERSP     = ERSP;
                study(end).ITC      = ITC;
                
            case 'both'
                study(end).Freqs    = frex;
                study(end).Data     = ERP;
                study(end).ERSP     = ERSP;
                study(end).ITC      = ITC;
        end
    else
        study(c_ind).ID(end+1,:) = {ID,ntrials};
        
        switch lower(type)
            case 'erp'
                study(c_ind).Data(:,:,end+1)     = ERP;
            case 'wavelet'
                study(c_ind).ERSP     = ERSP;
                study(c_ind).ITC      = ITC;
                
            case 'both'
                study(c_ind).Data(:,:,end+1)     = ERP;
                study(c_ind).ERSP(:,:,:,end+1)     = ERSP;
                study(c_ind).ITC(:,:,:,end+1)      = ITC;
        end
    end
    
    clear EEG ID cond ntrials data ERSP ITC c_ind frex
end

study = stucy(2:end);

%% save?

end