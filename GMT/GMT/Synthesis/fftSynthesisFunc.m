% function wav = fftSynthesisFunc(spec, par)
% Synthesize waveform from short-term Fourier coefficients.
% 
% INPUT:
%   spec  - matrix of complex Fourier coefficients; dimensions NFFT/2 x nFrames
%   par   - parameter struct or object
%
% FIELDS FOR PAR:
%   parent.nFFT - FFT length
%   parent.nHop - FFT hop size
%   parent.windowType - FFT analysis window type
%   synWindowType - synthesis window type, cf. generateWindow
%   synScaling - output scale factor ([] for auto-scale)
%                auto-scale: scale output such that the signal amplitude is
%                maintained through FFT analysis and subsequent re-synthesis
%                (for the chosen window types, frame overlap etc.) if no
%                further processing of the FFT coefficients occurs in between.
%   combineDcNy - DC and Nyquist components combined into single complex 1st bin? [0/1]
%                 If 1, then DC := Re(X1)+Im(X1), and NY := RE(X1)-Im(X1) 
%   compensateFftLength - Multiply FFT coefficients by nFft/2? [boolean] 
%
%
% OUTPUT:
%   wav - synthesized real-valued signal vector
%
% If synScaling = [], a scaling will be chosen automatically such that, 
% for no further processing between FFT and IFFT, the the envelope power of the signal
% prior to the FFT is approximately preserved in the output of the IFFT. 
% Otherwise, the output of the IFFT is multiplied by synScaling.
%
% See also FftSynthesisUnit, generateWindow
%
% Copyright (c) 2012-2020 Advanced Bionics. All rights reserved.
function wav = fftSynthesisFunc(spec, par)

% check for required parameters
requiredFields = {'synWindowType','synScaling'};
checkParamFields(par, requiredFields);

nFrames = size(spec,2);
nCoeff = size(spec,1);
strat = par.parent;
nFft = strat.nFft;
nHop = strat.nHop;
windowType = strat.windowType;
gOutUser = 10^(par.outputGain/20);

if par.compensateFftLength 
    spec = spec * (nFft/2);
end

FFTLENGTHS_ENV = 10;
if isempty(par.synScaling) % auto-scaling requested?
   gOutScaling = computeSynScaling(nFft, nHop, generateWindow(windowType,nFft), generateWindow(par.synWindowType,nFft), FFTLENGTHS_ENV); 
else
   gOutScaling = par.synScaling;
end

% generate synthesis window
synWin = generateWindow(par.synWindowType, nFft);
synWin = synWin(:); % ensure column vector

% Are DC and NF components coded into single complex FFT bin 1 ? -> unwind
if par.combineDcNy
    R = real(spec(1,:));
    I = imag(spec(1,:));
    spec(1,:) = R+I;
    spec(nFft/2+1,:) = R-I;
end

wav = real(istft(spec, synWin, nFft, nHop)) * gOutScaling * gOutUser ;

end

function out = istft(S, synWin, NFFT, hop)
lenBuff = NFFT;
nChan = size(S, 3);
nFrames = size(S,2);
nCoeff = size(S,1);

out = zeros(hop * (nFrames-1) + lenBuff,nChan);

S_i = zeros(NFFT,nChan);
revInd = ceil(NFFT/2):-1:2;
for iFrame = 0:nFrames-1
    idxWav = iFrame*hop+1:iFrame*hop+lenBuff;
    
    S_i(1:nCoeff,:) = S(:,iFrame+1,:);
    S_i(nCoeff+1:floor(NFFT/2)+1,:) = 0; % handles both nFft/2 and nFft/2+1 cases 
    S_i(floor(NFFT/2)+2:end,:) = conj(S_i(revInd,:));
    
    buff_i = ifft(S_i) .* repmat(synWin(:), [1 nChan]);
    
    out(idxWav,:) = out(idxWav,:) + buff_i;
end

end

