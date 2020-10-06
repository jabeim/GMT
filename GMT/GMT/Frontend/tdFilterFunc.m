% y = tdFilterFunc(par, x)
% Filter single-channel time-domain input signal x using one or multiple 
% filters defined by the numerator and denominator coefficients of a
% rational transfer function.
%
% INPUT:
%    x - 1 x nSamp input time signal 
%
% FIELDS FOR PAR:
%    coeffNum   - nFilters x nCoeff matrix of numerator coefficients of the z-transfer function;
%                 (each row constitutes one channel)
%    coeffDenom - nFilters x nCoeff matrix of denominator coefficients 
%                 (each row constitutes one channel)
% OUTPUT:
%    y - nFilters x nSamp filtered output signals
%
% See also: filter.m
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.

function Y = tdFilterFunc(par, x)
    
   nCh = size(par.coeffNum,1);
    
   Y = zeros(nCh, length(x));

    for iCh = 1:nCh
        Y(iCh,:) = filter(par.coeffNum(iCh,:), par.coeffDenom(iCh,:), x);
    end

end

