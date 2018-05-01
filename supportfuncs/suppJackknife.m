% This function preforms a jackknife opperation on ERP data and is called
% internaly by measument supplumetary functions.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
14-04-2018  Function re-write for better preformance
19-12-2016  Add Jackknife adjustment as per: Smulders, F. T. (2010) 
25-11-2016  New function (written in MATLAB R2015a)
%}

function jmean = suppJackknife(mode,data,dim)
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