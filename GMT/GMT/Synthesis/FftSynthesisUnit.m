% FftSynthesisUnit < ProcUnit
% Weighted overlap-add audio synthesis from complex short-term Fourier
% coefficients.
%
% FftSynthesisUnit properties:
%   synWindowType - synthesis window type (cf. generateWindowFunc) {'blackhann'}
%   synScaling - linear output scaling factor, or: 
%                  'autoAlways'/[] - do auto-scaling each time run is called 
%                  'autoOnce' - do auto-scaling at creation time; to re-do
%                        scaling at a later time, use "updateScalingFactor"
%                        method
%                Auto-scaling computes the factor based on the analysis and
%                synthesis FFT and window parameters, such that the signal
%                power is maintained through FFT analysis and subsequent 
%                re-synthesis (if no further processing occurs in between)
%   combineDcNy - DC and Nyquist components combined into single complex 1st bin? [boolean] [false]
%                 If 1, then DC := Re(X1)+Im(X1), and NY := RE(X1)-Im(X1)
%   compensateFftLength - Multiply FFT coefficients by nFft/2? [boolean] [false]
%
% FftSynthesisUnit methods:
%   FftSynthesisUnit(parent, ID, synWin, scaling) - constructor
%   run - execute processing
%
% Input ports:
%   #1 - complex short-term Fourier coefficient matrix, (NFft/2+1) x nFrames
%
% Output ports:
%   #1 - synthesized real-valued signal vector
%
% See also fftSynthesisFunc, generateWindow, FFTFilterBankUnit
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.

classdef FftSynthesisUnit < ProcUnit
    properties (SetAccess = immutable)
       synWindowType = 'blackhann'; % synthesis window type
       combineDcNy = false; % DC and Nyquist bins combined into single complex 1st bin? [false]
       compensateFftLength = false; % Multiply FFT coefficients by nFft/2? [boolean] [false]
       outputGain = 0; % gain applied to output waveform [dB] [0]
    end
        
    properties (SetAccess = private) % (SetObservable not necessary, because modified flag is always true)
        synScaling = []; % output scale factor ([] for auto-scale)        
    end
    
    methods
        function obj = FftSynthesisUnit(parent, ID, synWin, scaling, compensateFftLength, combineDcNy, outputGain)
        % obj = FftSynthesisUnit(parent, ID, synWin, scaling, combineDcNy)
            obj = obj@ProcUnit(parent, ID, 1, 1);
            
            obj.synWindowType = synWin;
                        
            if nargin > 3
                if isscalar(scaling)
                    obj.synScaling = scaling;
                elseif strcmpi(scaling, 'autoOnce')
                    obj.synScaling = obj.updateSynScaling();
                elseif strcmpi(scaling, 'autoAlways') || isempty(scaling)
                    obj.synScaling = [];
                else
                    error('Unknown scaling type %s', scaling);
                end
            end
          
            if nargin > 4
                obj.compensateFftLength = compensateFftLength;
            end
            
            if nargin > 5 && ~isempty(combineDcNy)
                obj.combineDcNy = combineDcNy;
            end
            
            if nargin > 6
                obj.outputGain = outputGain;
            end
            
        end
        
        function run(obj)
            spec = obj.getInput(1);
            wav  = fftSynthesisFunc(spec, obj);       
            obj.setOutput(1, wav);
        end
                
        % factor = updateSynScaling(obj)
        % (Re-)Estimate WOLA scaling factor based on current unit/strategy
        % parameters, and set the unit's synScaling property accordingly.
        function factor = updateSynScaling(obj)

            nFft = obj.parent.nFft;
            nHop = obj.parent.nHop;
            anaWin = generateWindow(obj.parent.windowType,nFft);
            synWin = generateWindow(obj.synWindowType,nFft);
            
            factor = computeSynScaling(nFft, nHop, synWin, anaWin, 10);
            
            obj.synScaling = factor;
        end
        
        function set.synScaling(obj, scaling)
            obj.synScaling = scaling;
        end
        
    end
    
end