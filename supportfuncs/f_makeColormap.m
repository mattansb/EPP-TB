% PURPOSE:  build your own custom colormaps using main colors rgbcmyk
%
%
% FORMAT
% ------
% cmap = f_makeColormap(colors)
%
%
% INPUTS
% ------
% colors    - string (char) of color codes - any sequence of 'rgbcmywk'
%             representing different colors (such as 'b' for blue) is
%             acceptable.
%
% EXAMPLE
% -------
%   figure; contourf(peaks);
%   colorbar; colormap(f_makeColormap('wkrygcbmmmmmmm'))
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
13-11-2018  Complete re-write
01-05-2018  New function (written in MATLAB R2017a)
%}

function cmap = f_makeColormap(colors)

%% Setup colors

d = [1 1 1;...
    1 0 0;...
    0 1 0;...
    0 0 1;...
    0 1 1;...
    1 0 1;...
    1 1 0;...
    0 0 0];

color_table = array2table(d,...
    'VariableNames',{'r','g','b'},...
    'RowNames',{'w','r','g','b','c','m','y','k'});

%% Get colors
colors  = arrayfun(@(x) {x},colors);
ncolors = length(colors);

x       = linspace(1,255,ncolors);
xq      = 1:255;
RBGs    = color_table{colors,{'r' 'g' 'b'}};
RBGsq   = interp1(x,RBGs,xq);

cmap = sqrt(RBGsq); % for smoothing: https://youtu.be/LKnqECcg6Gw

end