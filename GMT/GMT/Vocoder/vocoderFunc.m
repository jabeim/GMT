function [audioOut,audioFs] = vocoderFunc(par,electrodogram)
% [audioOut,audioFs] = vocoderFunc(par,electrodogram)
%
% The vocoder function transforms a matrix of electrode pulses into an
% acoustic wavform, optionally saving the data to disk as a file. The
% vocoder is intended to be used to perceptualize the changes of your
% signal processing. 
% It is not intended to be optimized by contestants.
% Any changes made to the vocoder will not persist during the judging phase
% of the contest.
%
% DO NOT EDIT ANY PARAMETERS THAT HAVE BEEN HARD CODED BY DEFAULT!!!
%
% INPUTS:
%   electrodogram: a matrix of cochlear implant current pulses by channel
%   over time (nChan x nSamp). Current amplitudes should be specified in
%   microamps [uA].
%
% FIELDS FOR PAR:
%   audioFs - the desired audio output sampling frequency. This is not
%       exact and is rounded to match internal vocoder processing requirements,48000 [Hz]
%   saveAudioOutput - whether or not to save the resulting output to a file, false [bool]
%   audioOutputFile - used as the name of the file if saveAudioOutput istrue, '', [string]
%
%
% OUTPUTS
%   audioOut - the floating point audio data resulting from the vocoder signal processing
%   audioFs - the calulated audio sampling frequency useful for playback.
%       This can be slightly different than specified in par.audioFs


%% Initialization
% set up default values for essential parameters in PAR
defaultParNames = {'audioFs','saveAudioOutput','audioOutputFile'};
defaultParValues = {48000,true,''};


%% DO NOT EDIT
% These values are set to match processing for contest submission. Changing
% them will result in audio that will not match what is heard during the
% judging phase

% Initialize contest-compatible vocoder settings
par.captFs = 55556;
par.nCarriers = 20;
par.elecFreqs = [];
par.spread = [];
par.neuralLocsOct = [];
par.nNeuralLocs = 300;
par.MCLmuA = [];
par.TmuA = [];
par.tAvg = .005;
par.tauEnvMS = 10;
par.nl = 5;



% scan optional arguments and assign values based on the defaults
for i = 1:length(defaultParNames)
    if ~isfield(par,defaultParNames{i})
        par.(defaultParNames{i}) = defaultParValues{i};
    end
end

%% scale and preprocess electrodogram data
% scale2MuA = 500/par.resistorVal;
scale2MuA = 1;

switch class(electrodogram)
    case 'char'
        electrodeAmp = readmatrix(electrodogram);
        if length(size(electrodeAmp)) == 2 && any(size(electrodeAmp) == 16)
            if size(electrodeAmp,1) == 16
            else
                electrodeAmp = electrodeAmp';
            end
        else
            error(['Electrodogram has improper size. Detected size: [' num2str(size(electrodeAmp,1)) ',' num2str(size(electrodeAmp,2)) '], the first dimension should be length 16.'])
        end
    case 'double'
%         electrodeAmp = electrodogram(:,2:end);   % the first column of data is actually channel similarity to the template.
        electrodeAmp = electrodogram;
    otherwise
        error('wrong electrodogram input')
end

nElec = size(electrodeAmp,1);
elData = electrodeAmp*scale2MuA;
captTs = 1/par.captFs;


%% compute electrode locations in terms of frequency
if isempty(par.elecFreqs)
    elecFreqs = logspace(log10(381.5),log10(5046.4),nElec);
else
    if nElec ~= length(par.elecFreqs)
        error('# of electrode frequencys does not match recorded data!')
    else
        elecFreqs = par.elecFreqs;
    end
end
   
%% electric field spread curve
if isempty(par.spread)
    elecPlacement = ones(1,nElec);
    load 'spread.mat'
    spread.fOct = fOct;
    spread.amp = voltage;
else
    elecPlacement = spread.elecPlacement;
    for i = 1:length(spread.curve)
        spread(i).fOct = spread.curve(i).fOct;
        spread(i).amp = spread.curve(i).amp;
    end
end
%% Octave location of neural populations
if isempty(par.neuralLocsOct)
    neuralLocsOct = [log2(linspace(150,850,40)) linspace(log2(870),log2(8000),260)];
else
    neuralLocsOct = par.neuralLocsOct;
end

if isempty(par.nNeuralLocs)
    nNeuralLocs = 300;
else
    nNeuralLocs = par.nNeuralLocs;
end

neuralLocsOct = interp1(1:length(neuralLocsOct),neuralLocsOct,linspace(1,length(neuralLocsOct),nNeuralLocs));

%% tauEnvMS
if isempty(par.tauEnvMS)
    tauEnvMS = 10;
else
    tauEnvMS = par.tauEnvMS;
end

taus = tauEnvMS/1000;
alpha = exp(-1/(taus*par.captFs));

%% MCL and T Levels in MuA
if isempty(par.MCLmuA)
    MCLmuA = 500*ones(1,nElec)*1.2;
else
    if length(par.MCLmuA == nElec)
        MCLmuA = par.MCLmuA*1.2;
    elseif length(par.MCLmuA == 1)
        MCLmuA = par.MCLmuA*ones(1,nElec)*1.2;
    else
        error('wrong number of MCL levels')
    end
end

if isempty(par.TmuA)
    TmuA = 50*ones(1,nElec);
else
    if length(par.TmuA == nElec)
        TmuA = par.TmuA;
    elseif length(par.TmuA == 1)
        TmuA = par.TmuA*ones(1,nElec);
    else
        error('wrong number of T levels')
    end
end
%% time frame to average neural activity
% If too high: smeared, if too low: low frequencies not reconstructed by
% window

if isempty(par.tAvg)
    tAvg = ceil(.005/captTs)*captTs;
else
    tAvg = ceil(par.tAvg/captTs)*captTs;
end
mAvg = round(tAvg/captTs);
blkSize = mAvg;

%% audio sample frequency of output
if isempty(par.audioFs)
    audioFs = fix(ceil(tAvg*44100)/tAvg);
else
    audioFs = fix(ceil(tAvg*par.audioFs)/tAvg);
end

audioTs = 1/audioFs;
nAvg = round(tAvg/audioTs);
tWin = 2*tAvg;
nFFT = round(tWin/audioTs);

%charge 2 EF matrix
charge2EF = zeros(nNeuralLocs,nElec);
elecFreqOct = log2(elecFreqs);

for iEl = 1:nElec  
    steerVec = interp1(...
        spread(elecPlacement(iEl)).fOct + elecFreqOct(iEl),...
        spread(elecPlacement(iEl)).amp,...
        neuralLocsOct,...
        'linear','extrap');
    steerVec(steerVec<0)=0;
    charge2EF(:,iEl) = steerVec(:);
end
%% matrix to map neural activity to FFT bin energies ||||||||||||||||||||||
MneurToBin=NeurToBinMatrix(neuralLocsOct,nFFT,audioFs);


%% other auxiliary variables
phs=2*pi*rand(floor(nFFT/2),1);%|||||||||||||||||||||||||||||||||||||||||||

audioPwr=zeros(nNeuralLocs,blkSize+1);

M=interp1(elecFreqOct,MCLmuA,neuralLocsOct);M=M(:);
M(neuralLocsOct<elecFreqOct(1))=MCLmuA(1);
M(neuralLocsOct>elecFreqOct(nElec))=MCLmuA(nElec);

T=interp1(elecFreqOct,TmuA,neuralLocsOct);T=T(:);
T(neuralLocsOct<elecFreqOct(1))=TmuA(1);
T(neuralLocsOct>elecFreqOct(nElec))=TmuA(nElec);

normRamp=bsxfun(@times,charge2EF,1./(M-T));
normOffset=T./(M-T);

elData(elData<0)=0;%if abs is used, ireelevant electrodes add constructively

% generate tone complex
nBlocks = (nFFT/2*floor(size(elData,2)/blkSize+1))-1;
tones = zeros(par.nCarriers,nBlocks);
toneFreqs = generate_cfs(20,20000,par.nCarriers);
t =(1:nBlocks)/audioFs;

for toneNum = 1:par.nCarriers
    tones(toneNum,:) = sin(2*pi*toneFreqs(toneNum)*t+phs(toneNum));
end

interpSpect = zeros(par.nCarriers,floor(size(elData,2)/blkSize));
fftFreqs = [1:floor(nFFT/2)]*audioFs/nFFT;

%% Loop through frames of electrode data, convert to electric field, calculate neural spectrum
for blkNumber=1:floor(size(elData,2)/blkSize)
    %% charge to patient-normalized electric field
    timeIdx=((blkNumber-1)*blkSize+1) : (blkNumber*blkSize);
    EF=max(0,bsxfun(@plus,normRamp*elData(:,timeIdx),-normOffset));
    
    EF = EF / 0.4 * 0.5;
    
    % normalized EF to activity
    nl = par.nl;
    % activity = max(0, min(exp(-nl+nl*EF),1)    - exp(-nl))/ (1 - exp(-nl));  
    activity = max(0, min(exp(nl*EF), exp(nl)) - 1 ) / (exp(nl) - 1);

    
    for k=1:blkSize
        audioPwr(:,k+1)=max(audioPwr(:,k)*alpha + activity(:,k)*(1-alpha),activity(:,k));
    end
    audioPwr(:,1)=audioPwr(:,end);    
    
    % averaged energy
    energy=sum(audioPwr,2)/mAvg;
    
    % spectrum of audio from average window
%     spect=(MneurToBin*energy).*exp(1i*phs);
%     
%     % reduce spectrum to tone frequency components
%     toneMags   = interp1(fftFreqs,abs(spect),toneFreqs,'linear','extrap');
%     tonePhases = interp1(fftFreqs,angle(spect),toneFreqs,'linear','extrap');
%     
%     interpSpect(:,blkNumber) = toneMags.*exp(1j*tonePhases);
    
    toneMags   = interp1(fftFreqs, (MneurToBin*energy), toneFreqs,'linear','extrap');
%     toneMags = (MneurToBin*energy);
    
    interpSpect(:,blkNumber) = toneMags;
end

    %interpolated spectral envelope scaling
    specVec = (1:blkNumber)*nFFT/2;
    newTimeVec = 0:(nBlocks-(nFFT/2));
%     interpSpect2 = zeros(length(toneFreqs),length(newTimeVec));
    modTones = zeros(length(toneFreqs),length(newTimeVec));
    
    %interpolate spectral envelope from frames back to td
    for freq = 1: length(toneFreqs)
%         tEnvMag = interp1(specVec,abs(interpSpect(freq,:)),newTimeVec,'linear','extrap');
%         tEnvPhs = interp1(specVec,angle(interpSpect(freq,:)),newTimeVec,'linear','extrap');
%         interpSpect2(freq,:) = tEnvMag.*exp(1j*tEnvPhs);
%         modTones(freq,:) = tones(freq,1:end-(nFFT/2-1)).*abs(interpSpect2(freq,:));
%         
        tEnvMag = interp1(specVec,interpSpect(freq,:),newTimeVec,'linear','extrap');
        modTones(freq,:) = tones(freq,1:end-(nFFT/2-1)) .* tEnvMag;
    end
    
    audioOut = sum(modTones,1);
    
    % fixed output level to -25dB so that differences in overall loudness
    % will not impact judging
    audioOut = audioOut/rms(audioOut)*10^(-25/20);
    
    if isfield(par,'saveAudioOutput') && par.saveAudioOutput == true
        if isfield(par,'audioOutputFile') && ~isempty(par.audioOutputFile)
            audiowrite(['Output' filesep par.audioOutputFile],audioOut,audioFs)
        else
            timestr = datestr(now,'yyyymmdd_HHMMSS');
            outFn = fullfile(['Output' filesep 'vocoderOutput_' timestr '.wav']);
            audiowrite(outFn,audioOut,audioFs)
        end
    end
    
end

