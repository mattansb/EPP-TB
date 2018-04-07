% This function preforms a jackknife opperation on ERP data and is called
% internaly by measument supplumetary functions.
%
% Author: Mattan S. Ben Shachar, BGU, Israel
%
% See also epp_getamplitude & epp_getamplitude

%{
Change log:
-----------
19-12-2016  Add Jackknife adjustment as per: Smulders, F. T. (2010) 
25-11-2016  New function (written in MATLAB R2015a)
%}

function dataOut = suppJackknife(d,dataIn)
switch d
    case 'in'
        [E, T, N] = size(dataIn);

        dataOut = zeros(E, T, N);
        nList = 1:N;

        for s = 1:N
            dataOut(:,:,s) = mean(dataIn(:,:,nList~=s),3);
        end
    case 'out'
        for ob = 1:size(dataIn,2) % each electrode
            N = length(dataIn(:,ob));
            J = nanmean(dataIn(:,ob));

            dataOut(:,ob) = N*J-(N-1)*dataIn(:,ob);
        end
end

end