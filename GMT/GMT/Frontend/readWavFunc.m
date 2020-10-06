% signalIn = readWavFunc(par)
%
% Read wav data from file. 
%
% INPUT:
%   par  - parameter object / struct
%
% FIELDS FOR PAR:
%   parent.fs - desireed output sampling frequency (Hz)
%   wavFile - name of wav-file (used of readWavFunc is called with a single arguement)
%   tStartEnd - 2-el vector defining start and end time of section to be read; sec
%   iChannel - index of channel to be returned [integer]
%
% OUTPUT:
%   signalIn - samples of wav-file
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.
function signalIn = readWavFunc(par)

name = par.parent.wavFile;
stratFs = par.parent.fs;

v = version; % get Matlab version
if sscanf(v,'%d') < 8  % before Matlab 2012b -> wavread
    [signalIn, srcFs] = wavread(name);    
else
    
    [path body ext] = fileparts(name);
    if isempty(ext)
        name = [name, '.wav'];
    end
    
    [signalIn, srcFs] = audioread(name,'native');  % 2012b or later -> audioread
    
    switch class(signalIn)
        case 'uint8'
            error('8 bit uint wav format not supported')
        case 'int16'
            bits = 16;
            maxBit = 2^(bits-1);
        case 'int32'
            bits = 32;
            maxBit = 2^(bits-1);
        case 'single'
            maxBit = 0;
        case 'double'
            maxBit = 0;
        otherwise
            error(['Data type ' class(signalIn) ' not supported.'])
    end
    signalIn = cast(signalIn,'double')/(maxBit+1);
    
end

signalIn = signalIn(:,par.iChannel);

if ~isempty(par.tStartEnd)
    iStartEnd = round(par.tStartEnd(:)*srcFs) + [1; 0];
    signalIn = signalIn(iStartEnd(1) : iStartEnd(2));
end

if srcFs ~= stratFs
    signalIn = resample(signalIn, stratFs, srcFs); 
end

