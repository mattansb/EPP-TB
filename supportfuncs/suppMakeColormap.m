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

%% Validate and initiate
p = inputParser;
    addOptional(p,'colors','wrgbcmyk',@ischar);
parse(p, colors); % validate

%% Init
ncolors = length(colors);
nbins   = ncolors-1;
bin_len = round((255/nbins)+2);

%% Colors

vec = cell([nbins 1]);

for i = 1:nbins
    R = linspace(...
        translateColor(colors(i),  'r'),...
        translateColor(colors(i+1),'r'),...
        bin_len...
        );
    
    G = linspace(...
        translateColor(colors(i),  'g'),...
        translateColor(colors(i+1),'g'),...
        bin_len...
        );
    
    B = linspace(...
        translateColor(colors(i),  'b'),...
        translateColor(colors(i+1),'b'),...
        bin_len...
        );
    
    vec{i} = [R; G; B]';
    
    if ~(i==1 || i==nbins)
        vec{i} = vec{i}(2:end-1,:);
    end
end


cmap = cat(1,vec{:});

end %end of buildcmap

function d = translateColor(char,place)
switch char
    case 'w'
        d = [1 1 1];
    case 'r'
        d = [1 0 0];
    case 'g'
        d = [0 1 0];
    case 'b'
        d = [0 0 1];
    case 'c'
        d = [0 1 1];
    case 'm'
        d = [1 0 1];
    case 'y'
        d = [1 1 0];
    case 'k'
        d = [0 0 0];
end

switch place
    case 'r'
        d = d(1);
    case 'g'
        d = d(2);
    case 'b'
        d = d(3);
end
end