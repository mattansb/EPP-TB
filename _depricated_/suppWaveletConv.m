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
% cycles_range  - [min max] | [num] A range or a constant of the number of
%                 cycels for each wavelet (also marked as the ratio between
%                 f0/sigma(f), or the n in n/(2*pi*f).)
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
21-05-2018  Fix bug when baseline correcting.
11-05-2018  Added support for constant f/sigma(f) / number of cycles in a
            wavelet.
10-05-2018  Added new baseline correction methods
13-04-2018  Added suuport for no time cuts
01-03-2018  Added support for downsampling.
28-02-2018  absolute ITC values are not computed,to allow for the
            combination of conditions.
13-03-2017  New function (written in MATLAB R2015a)
%}

function [new_power,itpc,frex,cut_times] = suppWaveletConv(EEG,freq_range,num_frex,cycles_range,baselinetime,cut_times,varargin)

[new_power,itpc,frex,cut_times] = f_WaveletConv(EEG,freq_range,num_frex,cycles_range,baselinetime,cut_times,varargin{:});
warning('suppWaveletConv is no longer supported. Use f_WaveletConv instead.')

end
