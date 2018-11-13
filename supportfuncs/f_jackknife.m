% PURPOSE:  Jackknife transform and back-trasnform a matrix.
%
%
% FORMAT
% ------
% jmean = f_jackknife(mode,data,dim)
%
%
% INPUTS
% ------
% mode          - ['in'|'out'] Whether data should be transformed or
%                 back-transformed.
% data          - A matrix of any size and any number of dimentions.
% dim           - The dimentins along which to preform the jackknife.
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
14-04-2018  Function re-write for better preformance
19-12-2016  Add Jackknife adjustment as per: Smulders, F. T. Y. (2010) (DOI: 10.1111/j.1469-8986.2009.00934.x)
25-11-2016  New function (written in MATLAB R2015a)
%}

function jmean = f_jackknife(mode,data,dim)
switch lower(mode)
    case 'in'
        % Get size
        mat_size    = size(data);
        r_size      = ones([1 length(mat_size)]);
        r_size(dim) = mat_size(dim);

        % Get mean
        M = repmat(mean(data,dim), r_size);

        % Get Jacked mean
        jmean = (M*r_size(dim)-data)/(r_size(dim)-1);
    case 'out'        
        % Get size
        mat_size    = size(data);
        r_size      = ones([1 length(mat_size)]);
        r_size(dim) = mat_size(dim);

        % Get mean
        J = repmat(mean(data,dim), r_size);
        
        % Get Un-Jacked mean
        jmean = J*r_size(dim)-data*(r_size(dim)-1);
end
end