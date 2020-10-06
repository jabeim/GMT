% FftFilterbankUnit < ProcUnit
% Compute FFT on buffered signal segments.
%
% FftFilterbankUnit properties
%   combineDcNy - Combine DC and Nyquist bins into single complex 1st bin? [boolean] [false]
%                 if true, then bin #1 = .5*(DC+NY + i(DC-NY)) and bin #NFFT/2+1 = 0 
%   compensateFftLength - Divide FFT coefficients by nFft/2? [boolean] [false] 
%   includeNyquistBin - return bin #nFft/2+1 in output? [boolean] [false]  
% 
% Input Ports:
%   #1  - buffered signal frames, NFFT x nFrames
%
% Output Ports:
%   #1  - FFT coefficient matrix, (NFFT/2) x nFrames (default) or (NFFT/2+1) x nFrames, depending on includeNyquistBin
%
% See also: fftFilterBankFunc.m
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.
classdef FftFilterbankUnit < ProcUnit
    
    properties (SetAccess = private)
        combineDcNy = false; % Combine DC and Nyquist bins into single complex 1st bin?
        compensateFftLength = false; % Divide FFT coefficients by nFft/2? [boolean]
        includeNyquistBin = false; % Return bin #nFft/2+1 in output? [boolean] 
    end
   
    methods
        function obj = FftFilterbankUnit(parent, ID, compensateFftLength, combineDcNy, includeNy)
        % obj = FftFilterbankUnit(parent, ID, compensateFftLength, combineDcNy)
            obj = obj@ProcUnit(parent, ID, 1, 1);
            if nargin > 2
                obj.compensateFftLength = compensateFftLength;
            end

            if nargin > 3
                obj.combineDcNy = combineDcNy;
            end
            if nargin > 4
                obj.includeNyquistBin = includeNy;
            end
        end
        function run(obj)
            buf = obj.getInput(1);
            X = fftFilterbankFunc(obj, buf);
            obj.setOutput(1, X);
        end
    end
end