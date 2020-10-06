% X = fftFilterbankFunc(par, buf)
% Compute FFT along each column of a matrix buffered signal segments.
%
% INPUT:
%   par - parameter object / struct
%   buf - nFft x nFrames matrix of windowed signal buffers
%
% FIELDS FOR PAR:
%   parent.nFft - FFT length [int > 0]
%   combineDcNy - Combine DC and Nyquist bins into single complex 1st bin? [boolean]
%                 if true, then bin #1 = .5*(DC+NY + i(DC-NY)) and bin #NFFT/2+1 = 0
%   compensateFftLength - Divide FFT coefficients by nFft/2? [boolean]
%   includeNyquistBin - return bin #nFft/2+1 in output? [boolean]
%
% OUTPUT:
%   X - FFT coefficient matrix, (NFFT/2) x nFrames or (NFFT/2+1) x nFrames depending on includeNyquistBin
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.
function X = fftFilterbankFunc(par, buf)

nFft = par.parent.nFft;

X = fft(buf, nFft);

if par.combineDcNy
    NF = X(nFft/2+1, :, :);
    DC = X(1, :, :);
    
    X(1, :, :) = (real(DC) + real(NF)) + 1j*(real(DC) - real(NF));    
    X(nFft/2+1, :, :) = 0;
end

if par.includeNyquistBin
    X = X(1:nFft/2+1, :, :); % positive frequencies incl. nyquist
else
    X = X(1:nFft/2, :, :); % positive frequencies excl. nyquist
end    

if par.compensateFftLength
    X = X ./ (nFft/2);
end