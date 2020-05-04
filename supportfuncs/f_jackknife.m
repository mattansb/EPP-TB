% PURPOSE:  Jackknife transform and back-trasnform a matrix.
%
%
% FORMAT
% ------
% [results, WM_tot] = f_jackknife('in', data,dim,W)
%          results  = f_jackknife('out',data,dim,W,WM)
%
%
% INPUTS
% ------
% data      - A matrix of any size and any number of dimentions, to compute
%             jackknifed means from, or to back trasnsform according to
%             Smulders, F. T. Y. (2010) (DOI:
%             10.1111/j.1469-8986.2009.00934.x)
% dim       - The dimentins along which to preform the jackknife / back
%             transform.  
% W         - (optional) A vector of weights for computing the jackknife
%             means / back transforming the values. If left blank
%             un-weighted jackknife means / back transformed values.
% WM        - (optional) A values to be used as the mean value when back
%             trasforming jackknifed values (See Sluth & Meiran,
%             https://doi.org/10.1101/403766). If left blank, the inversing
%             is done around the mean jackknifed value.
%
%
% OUTPUT(s)
% ---------
% results   - A matrix the same size as data, of either jackknifed means,
%             or inversed jackknifed values.
% WM_tot    - For computing the jackknifed means, the total (weighted /
%             unweighted) mean. This can be used to computed the desired
%             value, which can then be used as the WM input when back
%             transforming the jackknife values.
%
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%{
Change log:
-----------
01-12-2018  Added support of weighted jackknifing
            Added support for custom central value See Sluth & Meiran
            https://doi.org/10.1101/403766.
14-04-2018  Function re-write for better preformance
19-12-2016  Add Jackknife adjustment as per: Smulders (2010),
            https://doi.org/10.1111/j.1469-8986.2009.00934.x 
25-11-2016  New function (written in MATLAB R2015a)
%}

function [Xout, WMout] = f_jackknife(mode,Xin,dim,W,WM)
% for backward compatability
if exist('W','var')~=1,     W   = [];   end
if exist('WM','var')~=1,    WM  = [];   end

% Set some defaults
if isempty(dim),    dim = find(size(Xin)>1, 'first');   end
if isempty(W),      W   = 1;                            end

% Get weights
W = resize_weights(W,Xin,dim);

switch lower(mode)
    case 'in'
        % Get (weighted) mean
        WMout   = sum(W.*Xin,dim)./sum(W,dim);
        WM      = resize_mean(WMout,Xin,dim);

        % Get jackknifed (weighted) means
        Xout = (WM.*sum(W,dim)-Xin.*W)./(sum(W,dim) - W);
    case 'out'
        % Get (weighted) mean
        if isempty(WM), WM = mean(Xin,dim); end
        WM      = resize_mean(WM,Xin,dim);

        % Get un-jackknifed (weighted) values
        Xout = (WM.*sum(W)-Xin.*(sum(W) - W))./W;
end
end

function W = resize_weights(W,X,dim)
if length(W)==1
    size_W      = ones([1 length(size(X))]);
    size_W(dim) = size(X,dim);
    W           = repmat(W,size_W);
elseif size(W,dim) ~= size(X,dim)
    if length(W)==size(X,dim)
        order_W         = 1:length(size(X));
        order_W(dim)    = 1;
        order_W(1)      = dim;
        W               = permute(W,order_W);
    else
        error('dimentional length of W must match specified dimention of X')
    end
end
end

function WM = resize_mean(WM,X,dim)
size_WM         = ones([1 length(size(X))]);
size_WM(dim)    = size(X,dim);
WM              = repmat(WM,size_WM);
end