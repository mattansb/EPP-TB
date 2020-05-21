ERPS = epp_loadeeglab(true, 'erp', true, 'combine', true);

t = f_timewindowIndex([-200 600], 'times', ERPS(1).timeLine);
for c = 1:length(ERPS)
    ERPS(c).timeLine = ERPS(c).timeLine(t);
    ERPS(c).Data = ERPS(c).Data(:,t,:);
end

WAV = epp_loadeeglab(true, 'wavelet', true, 'combine', true);

load('doc\sample data\chanlocs.mat')


%% butterfly

epp_plotbutterfly(ERPS, {ERPS.Condition}, [4 5 6])

epp_plotbutterfly(WAV, {ERPS.Condition}, [4 5 6], ...
    'freqs',[3 5; 5 10; 10 20], 'type', 'ersp')

epp_plotbutterfly(WAV, {ERPS.Condition}, [4 5 6], ...
    'freqs',[3 5; 5 10; 10 20], 'type', 'itc')

%% trace

epp_plotbutterfly(ERPS, {ERPS.Condition}, [4 5 6], 'trace', true)

epp_plotbutterfly(WAV, {ERPS.Condition}, [4 5 6], 'trace', true, ...
    'freqs',[3 5; 5 10; 10 20], 'type', 'ersp')

epp_plotbutterfly(WAV, {ERPS.Condition}, [4 5 6], 'trace', true, ...
    'freqs',[3 5; 5 10; 10 20], 'type', 'itc')

%% chanels

epp_plotchannels(ERPS, {ERPS.Condition}, [], 'chanlocs', chanlocs)

epp_plotchannels(WAV, {ERPS.Condition}, [], 'chanlocs', chanlocs,...
    'freqs',[3 5; 5 10; 10 20], 'type', 'ersp')

epp_plotchannels(WAV, {ERPS.Condition}, [], 'chanlocs', chanlocs,...
    'freqs',[3 5; 5 10; 10 20], 'type', 'itc')


%% Grands

epp_plotgrands(ERPS, {ERPS.Condition}, [4 5 6], 'errorType', 'SE')

epp_plotgrands(WAV, {ERPS.Condition}, [4 5 6], 'errorType', 'SE', ...
    'freqs',[3 5; 5 10; 10 20], 'type', 'ersp')

epp_plotgrands(WAV, {ERPS.Condition}, [4 5 6], 'errorType', 'SE',...
    'freqs',[3 5; 5 10; 10 20], 'type', 'itc')


%% TF

epp_plotTF(WAV, {ERPS.Condition}, [4 5 6])

%% topo

epp_plottopo(ERPS, chanlocs, {ERPS.Condition}, [0 100 200])

epp_plottopo(WAV, chanlocs, {ERPS.Condition}, [0 100 200],...
    'freqs',[3 5; 5 10; 10 20], 'type', 'itc')

epp_plottopo(WAV, chanlocs, {ERPS.Condition}, [0 100 200],...
    'freqs',[3 5; 5 10; 10 20], 'type', 'ersp')




