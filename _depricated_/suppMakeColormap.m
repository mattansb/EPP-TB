% PURPOSE:  build your own custom colormaps using main colors rgbcmyk
%
%
% FORMAT
% ------
% cmap = suppMakeColormap(colors)
%
%
% INPUTS
% ------
% colors    - string (char) of color codes - any sequence of rgbcmywk
%             representing different colors (such as 'b' for blue) is
%             acceptable.
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
01-05-2018  New function (written in MATLAB R2017a)
%}

function cmap = suppMakeColormap(colors)

cmap = f_makeColormap(colors);
warning('suppMakeColormap is no longer supported. Use f_makeColormap instead.')

end