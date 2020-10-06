function b = winBufFunc(par, signalIn)
% b = winBufFunc(par, signalIn)
% INPUT:
%   signalIn - samples of input signal obtained from e.g. wavread
%
% FIELDS FOR PAR:
%   parent.nFft - buffer size in samples (probably equal to FFT size)
%   parent.nHop - hop size in samples
%   parent.window - samples of window function
%   bufOpts - initial buffer state prior to signal onset. 'nodelay' start
%            buffering with first full input frame; [] for leading zeros; 
%            vector of length (nFft-nHop) to define arbitrary state.
%
% OUTPUT:
%   buf - buffers, one signal-frame per column
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.

strat = par.parent;

[M, N] = size(signalIn);
if N>M
%    warning('winBufFun: input signal wider than long (%dx%d). Transposing.', M, N);
    signalIn = signalIn.';
    N=M;
end

b = buffer(signalIn(:, 1), strat.nFft, strat.nFft-strat.nHop, par.bufOpt);
b = bsxfun(@times, b, strat.window);
if N>1
    temp = b;
    b = zeros(size(b, 1), size(b, 2), N);
    b(:, :, 1) = temp;
    
    for n=2:N
        b(:, :, n) = buffer(signalIn(:, n), strat.nFft, strat.nFft-strat.nHop, par.bufOpt);
        b(:, :, n) = bsxfun(@times, b(:, :, n), strat.window);
    end
end
