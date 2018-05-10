% PURPOSE:  compute ersp and itc.
%
%
% FORMAT
% ------
% [new_power,itpc,frex,cut_times] = suppWaveletConv3(EEG,freq_range,num_frex,cycles_range,baselinetime,cut_times,varargin)
%
%
% INPUTS
% ------
% EEG           - eeglab data struct (will compute for all channels)
% freq_range    - [min max] frequency (Hz) range to compute
% num_frex      - number of frequencies to compute (NOTE: will be]
%                 log-scaled).
% cycles_range  - [min max] range of cycels to use for range of
%                 frequencies.
% baselinetime  - [start end] baseline (in ms) for baseline correction
% cut_times     - [start end] range of time to save. if empty, all time
%                 point are returned.
% 
% Optional arguments
% ------------------
% 'baseline'    - one of the following baseline correction methods:
%                 'dB' (defult) - 10*log10(power/mean_baseline_power)
%                 'normalize'   - 100*(power - mean_baseline_power)/mean_baseline_power
%                 'standardize' - (power - mean_baseline_power)/std_baseline_power
%                 (see http://neuroimage.usc.edu/brainstorm/Tutorials/TimeFrequency#Morlet_wavelets)
% 'log'         - if true [defult], frquenceis are taken from logspace. If
%                 false, frequencies are taken from linspace.
% 'downsample'  - factor by witch to down sample (defult: 1 - no
%                 downsampling). See downsample().
% 'sound'       - ['big' (defult) |'beep'|'off'] sound to play when done.
% 
%
%
% Author: Mattan S. Ben Shachar & Rachel Rac, BGU, Israel

%{
Change log:
-----------
10-05-2018  Added new baseline correction methods
13-04-2018  Added suuport for no time cuts
01-03-2018  Added support for downsampling.
28-02-2018  absolute ITC values are not computed,to allow for the
            combination of conditions.
13-03-2017  New function (written in MATLAB R2015a)

TO DO
add option of '1/f compensation' by multiplying each power by the frequancy.
this step should happen before baseline correction (Baseline correction negates
this step - as the baseline correction cancels out this correction)
This is called  Spectral flattening
*** Would only be relevant for standardize method. All other methods it
mitkazez
%}

function [new_power,itpc,frex,cut_times] = suppWaveletConv3(EEG,freq_range,num_frex,cycles_range,baselinetime,cut_times,varargin)

tic
%% Validate
p = inputParser;
    addRequired(p,'EEG',@isstruct);
    addRequired(p,'freq_range',@(x) length(x)==2 && isnumeric(x));
    addRequired(p,'num_frex',@(x) length(x)==1 && isnumeric(x));
    addRequired(p,'cycles_range',@(x) length(x)==2 && isnumeric(x));
    addRequired(p,'baselinetime',@(x) length(x)==2 && isnumeric(x));
    addOptional(p,'cut_times',[],@(x) isempty(x) | (length(x)==2 && isnumeric(x)));
    addParameter(p,'log',true, @islogical);
    addParameter(p,'baseline','dB', @(x) any(strcmpi(x,{'dB','normalize','standardize'})));
    addParameter(p,'dB',nan, @islogical); % deprecated 
    addParameter(p,'downsample',1, @(x) floor(x) == x);
    addParameter(p,'sound','big', @(x) any(strcmpi(x,{'big','beep','off'})));
parse(p ,EEG,freq_range,num_frex,cycles_range,baselinetime,cut_times,varargin{:}); % validate


%% Define Parameters

nbchan          = size(EEG.data,1); % number of channels
min_freq        = freq_range(1);    % min freq to compute
max_freq        = freq_range(2);    % max freq to compute
num_frex        = num_frex;         % number of frex to compute
wavelet_cycles1 = cycles_range(1);  % for min freq
wavelet_cycles2 = cycles_range(2);  % for max freq

% Length of wavelet in seconds (p. 167-168)
time    = -2:1/EEG.srate:2;

% Frequencies (Hz) & number of cycles for each (p. 168, 196)
if p.Results.log % LOG-SCALED (recomended)
    frex    = logspace(log10(min_freq),log10(max_freq),num_frex);
    cyc     = logspace(log10(wavelet_cycles1),log10(wavelet_cycles2), num_frex)./(2*pi*frex); 
else % linear-scale
    frex    = linspace(min_freq,max_freq,num_frex);
    cyc     = linspace(wavelet_cycles1,wavelet_cycles2, num_frex)./(2*pi*frex);
end

%% Power & ITPC

% Definte convolution parameters:
n_wavelet            = length(time);
n_data               = EEG.pnts*EEG.trials;
n_convolution        = n_wavelet+n_data-1;
n_conv_pow2          = pow2(nextpow2(n_convolution));
half_of_wavelet_size = (n_wavelet-1)/2;


% initialize
eegfft      = zeros(nbchan,n_conv_pow2);
itpc        = zeros(nbchan,length(frex),EEG.pnts);
eegpower    = zeros(nbchan,num_frex,EEG.pnts); % elec X frequencies X time
new_power = zeros(nbchan,num_frex,EEG.pnts);


% get FFT of data
for ch = 1:nbchan % each channel
    eegfft(ch,:) = fft(reshape(EEG.data(ch,:,:),1,EEG.pnts*EEG.trials),n_conv_pow2);
end


% Loop through frequencies and compute synchronization
fprintf('Computing power and coherence... ');
pc100   = 0;
pp      = fprintf('%1$d%%',pc100);
for ch = 1:nbchan % each channel
    for fi = 1:num_frex % each freq        
        % Create Wavelet for current frequency
        wavelet = fft( sqrt(1/(cyc(fi)*sqrt(pi))) * exp(2*1i*pi*frex(fi).*time) .* exp(-time.^2./(2*(cyc(fi)^2))) , n_conv_pow2 );

        % Perform convolution
        eegconv = ifft(wavelet.*eegfft(ch,:));
        eegconv = eegconv(1:n_convolution);
        eegconv = eegconv(half_of_wavelet_size+1:end-half_of_wavelet_size);

        % Average power over trials 
        temppower           = mean(abs(reshape(eegconv,EEG.pnts,EEG.trials)).^2,2);
        eegpower(ch,fi,:)   = temppower;

        % extract ITPC
	% OR: don't abs so that conditions can be averaged. add abs into processing functions (plots, measures... etc)
        eegconv_ITPC    = reshape(eegconv,EEG.pnts,EEG.trials);
%         temp_ITPC       = abs(mean(exp(1i*angle(eegconv_ITPC)),2));
        temp_ITPC       = mean(exp(1i*angle(eegconv_ITPC)),2);
        itpc(ch,fi,:)   = temp_ITPC;
        
        tmp_pc = round(((ch-1)*num_frex+fi)*100/(num_frex*nbchan));
        if tmp_pc~=pc100 % use mywaitbar instead
            fprintf(repmat('\b',[1 pp]));
            pc100   = tmp_pc;
            pp      = fprintf('%1$d%%',pc100);
        end
    end
end
fprintf([repmat('\b',[1 pp]) 'done! '])

%% Baseline correction
if ~isnan(p.Results.dB)
    if double(p.Results.dB)==1
        warning('{''dB'',true} is no longer in use. Use {''baseline'',''dB''} instead.')
        p.Results.baseline = 'dB';
    elseif double(p.Results.dB)==0
        warning('{''dB'',false} is no longer in use. Using {''baseline'',''divide''}. The {''baseline'',''normalize''} should be used instead.')
        p.Results.baseline = 'divide';
    end    
end


[~,baselineidx(1)] = min(abs(EEG.times-baselinetime(1)));
[~,baselineidx(2)] = min(abs(EEG.times-baselinetime(2)));

switch p.Results.baseline
    case 'dB'
        fprintf('\nConverting power to dB... ');
        for ch = 1:nbchan % each channel
            temp_pow            = reshape(eegpower(ch,:,:),num_frex,EEG.pnts);
            baseline_power      = mean(temp_pow(:,baselineidx(1):baselineidx(2)),2);    % mean power in baseline
            new_power(ch,:,:)   = 10*log10(bsxfun(@rdivide,temp_pow,baseline_power));   % convert to dB
        end
    case 'normalize'
        fprintf('\nNormalizing power... ');
        for ch = 1:nbchan % each channel
            temp_pow            = reshape(eegpower(ch,:,:),num_frex,EEG.pnts);
            baseline_power      = mean(temp_pow(:,baselineidx(1):baselineidx(2)),2);    % mean power in baseline
            new_power(ch,:,:)   = 100*bsxfun(@minus,temp_pow,baseline_power);           % subtract mean
            new_power(ch,:,:)   = bsxfun(@rdivide,new_power(ch,:,:),baseline_power);    % devide by mean
        end
    case 'standardize'
        fprintf('\nStandardizing power... ');
        temp_pow            = reshape(eegpower(ch,:,:),num_frex,EEG.pnts);
            baseline_power      = mean(temp_pow(:,baselineidx(1):baselineidx(2)),2);    % mean power in baseline
            baseline_std        = std(temp_pow(:,baselineidx(1):baselineidx(2)),2);     % std power in baseline
            new_power(ch,:,:)   = 100*bsxfun(@minus,temp_pow,baseline_power);           % subtract mean
            new_power(ch,:,:)   = bsxfun(@rdivide,new_power(ch,:,:),baseline_std);      % devide by std
    case 'divide'
        fprintf('\nSemi-Normalizing power... ');
        for ch = 1:nbchan % each channel
            temp_pow            = reshape(eegpower(ch,:,:),num_frex,EEG.pnts);
            baseline_power      = mean(temp_pow(:,baselineidx(1):baselineidx(2)),2);    % mean power in baseline
            new_power(ch,:,:)   = bsxfun(@rdivide,temp_pow,baseline_power);             % devide by mean
        end
end
fprintf('done!');

%% Save selected time points
% Downsample the data, or remove buffer-zone segments.
fprintf('\nSaving... ');
if isempty(cut_times)
    cut_times = EEG.times;
else
    cut_times = EEG.times(EEG.times>=cut_times(1) & EEG.times<=cut_times(2));
end

% down sample?
if p.Results.downsample > 1
    cut_times = downsample(cut_times,p.Results.downsample);
end

% Save POWER
timesidx    = dsearchn(EEG.times',cut_times');
new_power   = new_power(:,:,timesidx);
itpc        = itpc(:,:,timesidx);

%% SOUND
switch p.Results.sound
    case 'big'
        load('handel.mat');
        sound(y(1000:18000))
    case 'beep'
        beep
end
fprintf('DONE. (elapsed time %d seconds)\n',round(toc))

end

% tic
% function [step, strLength] = mywaitbar(compl, total, step, nSteps, strLength)
% 
% progStrArray = '/-\|';
% tmp = floor(compl / total * nSteps);
% if tmp > step
%     fprintf(1, [repmat('\b', 1, strLength) '%s'], repmat('=', 1, tmp - step))
%     step = tmp;
%     ete = ceil(toc / step * (nSteps - step));
%     strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '%s %3d%%, ETE %02d:%02d'], progStrArray(mod(step - 1, 4) + 1), floor(step * 100 / nSteps), floor(ete / 60), mod(ete, 60));
% end
