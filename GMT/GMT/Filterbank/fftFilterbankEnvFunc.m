% function Env = fftFilterbankEnvFunc(X, mapFft2Ch)
% Compute the amplitude envelope at the output of a FFT filterbank for FFT
% coefficients X and filterbank channels defined by mapFft2Ch.
%
% INPUT:
%   X - FFT coefficients, size nFFT x nFrames
%   par - parameter struct or object
%
% FIELDS FOR PAR:
%   mapFft2Ch - EITHER:
%                 a 2-el. vector, specifying (1) the number of channels 
%                 for a F120 filterbank (3..15) and (2) extended low (0/1); 
%                 taylored for a fs=17400, nFFT=256
%               OR: a scalar, specifying the number of ch. for a F120
%                 filterbank without extended low channel
%               OR: an nCh x nFft mixing matrix (linear weights, >= 0)
%               OR: []  to use the strategy's startBin / nBinLims channel allocation
%   parent.startBin - first FFT bin in lowest filter  (needed if mapFft2Ch = [])
%   parent.nBinLims - number of FFT bins per channel  (needed if mapFft2Ch = [])
%
% OUTPUT:
%   Env - Filterbank envelopes, size nCh x nFrames
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.
function Env = fftFilterbankEnvFunc(X, par)
    
    strat = par.parent;
    startBin = strat.startBin;
    nBinLims = strat.nBinLims;
    
    if isempty(par.mapFft2Ch)
        mapFft2Ch = binLimsToMixingWeights(startBin, nBinLims, size(X,1));
    elseif isscalar(par.mapFft2Ch)
        mapFft2Ch = computeF120FilterbankWeights(par.mapFft2Ch, 0);
    elseif isvector(par.mapFft2Ch) % compute standard F120 FB weights
        mapFft2Ch = computeF120FilterbankWeights(par.mapFft2Ch(1),par.mapFft2Ch(2));
    else
        mapFft2Ch = par.mapFft2Ch;
    end

    assert(size(X,1) == size(mapFft2Ch,2), 'Matrix dimensions must match: size(X,1) == size(mapFft2Ch,2)');
    assert(all(mapFft2Ch(:) >= 0), 'Mixing weights must be non-negative')
    
    Env = sqrt(mapFft2Ch * (abs(X).^2));
end