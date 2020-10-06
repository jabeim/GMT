% WinBufUnit < ProcUnit
% Divide an input signal vector into overlapping buffers with a window
% function applied.
%
% WinBufUnit properties:
%  bufOpt - initial buffer state prior to signal onset. 'nodelay' start
%            buffering with first full input frame; [] for leading zeros; 
%            vector of length (nFft-nHop) to define arbitrary state.
%
% Input ports:
%   #1 - input signal vector (row or column). If matrix, assume it's a
%        multichannel signal where the longest dimension is time 
% Output ports:
%   #1 - matrix of buffered signal frames, NFFT x nFrames, or if
%        multichannel, a NFFT x nFrames x channels tensor.
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.

classdef WinBufUnit < ProcUnit

    properties (SetAccess = immutable)
        bufOpt = []; 
    end
    
    methods
        function obj = WinBufUnit(parent, ID, bufOpt)
            obj = obj@ProcUnit(parent, ID, 1, 1);
            if nargin > 2
                obj.bufOpt = bufOpt;
            end
        end
        
        function run(obj)
           signalIn = obj.getInput(1);
           buf = winBufFunc(obj, signalIn);
           obj.setOutput(1, buf);                      
        end
    end
end