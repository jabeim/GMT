% ReadWavUnit < ProcUnit
% 
% Read wav data from file. 
%
% Input Ports: none
%
% Output Ports:
%   #1  - column vector contained wav data
%
% ReadWavUnit Properties:
%  *wavFile - wav file name
%  *tStartEnd - 2-element vector specifying start and end time of the
%              section to be read [seconds]
%  *iChannel - index of channel to be returned [integer] [1]
%
% ReadWavUnit Methods:
%   ReadWavUnit - constructor
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.
classdef ReadWavUnit < ProcUnit
    
    properties (SetObservable, AbortSet)
        wavFile = '';  % wav file name
        tStartEnd = []; % 2-el vector defining start and end time of section to be read; sec
        iChannel = 1;
    end
    
    methods
        function obj = ReadWavUnit(parent, ID, wavFile, tStartEnd, iChannel)
        % obj = ReadWavUnit(parent, ID, wavFile, tStartEnd)
        % Constructor, generate new ReadWavUnit object
        
        obj = obj@ProcUnit(parent, ID, 0, 1);
               
            if (nargin > 2)
                obj.wavFile = wavFile;
            end
            
            if (nargin > 3)
                obj.tStartEnd = tStartEnd;
            end
            
            if (nargin > 4)
                obj.iChannel = iChannel;
            end
        end
        
        function run(obj) 
            signalIn = readWavFunc(obj);
            obj.setOutput(1, signalIn);
        end
    end
end