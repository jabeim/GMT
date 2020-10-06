% CurrentSteeringWeightsUnit < ProcUnit
% Compute current steering weights from nominal cochlear location of the 
% estimated spectral peak frequency per channel. Assumes that locations 
% for channel i lie within [i-1, i]. Steering weights can be distretized 
% and/or limited to a channel-specific sub-range. For each channel, the 
% resulting pairs of steering weights sum to 1 always. 
%
% CurrentSteeringWeightsUnit properties:
%   nDiscreteSteps - number of discretization steps [int >= 0] [9]; 0 -> no discretization
%   steeringRange  - range of steering between electodes [0..1] [1.0]
%                         - scalar range around 0.5 for all channels (within [0,1])
%                         - 1 x nChan vector of ranges (0-1) around 0.5 per channel
%                         - 2 x nChan matrix with (absolute) lower and upper steering 
%                           limits (0-1) per channel
%
% CurrentSteeringWeightsUnit methods:
%   CurrentSteeringWeightsUnit(parent, ID) - constructor
%   run() - execute processing
% 
% Input ports:
%   #1 - nChan x nFrames matrix of cochlear location, each row i limited to
%        values in [i-1, 1]
%
% Output ports:
%   #1 - (2*nChan) x nFrames matrix of current steering weights; weights 
%        for the lower and higher electrode of channel i are contained in
%        rows i and (i+nChan), resp.
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.

classdef CurrentSteeringWeightsUnit < ProcUnit
    properties (SetObservable)
        nDiscreteSteps = 9;   % nr. of discretization steps  [int >= 0] [9]; 0 -> no discretization
        steeringRange = 1.0;  % steering range between electrodes [0..1] [1.0]
    end
    
    methods
        function obj = CurrentSteeringWeightsUnit(parent, ID)
            obj = obj@ProcUnit(parent, ID, 1, 1);
        end
        
        function run(obj)
            loc = obj.getInput(1);
            weights = currentSteeringWeightsFunc(obj, loc);
            obj.setOutput(1, weights);            
        end
    end
end