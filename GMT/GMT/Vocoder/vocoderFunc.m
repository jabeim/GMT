function [audioOut,audioFs] = vocoderFunc(par,electrodogram)
% scan input parameters adding necessary default values
defaultParNames = {'captFs','nCarriers','elecFreqs','spread','neuralLocsOct','nNeuralLocs','MCLmuA','TmuA','tAvg','audioFs','tauEnvMS','nl','resistorVal'};
defaultParValues = {200000,20,[],[],[],300,[],[],.005,48000,10,8,2};

for i = 1:length(defaultParNames)
    if ~isfield(par,defaultParNames{i})
        par.(defaultParNames{i}) = defaultParValues{i};
    end
end

%% scale and preprocess electrodogram data
scale2MuA = 500/par.resistorVal;
electrodeAmp = electrodogram;
nElec = size(electrodogram,1);
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
    audioFs = ceil(tAvg*44100)/tAvg;
else
    audioFs = ceil(tAvg*par.audioFs)/tAvg;
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
    
    
    % normalized EF to activity
    nl = par.nl;
    activity=max(0,min(exp(-nl+nl*EF),1)-exp(-nl))/(1-exp(-nl));  
    
    for k=1:blkSize
        audioPwr(:,k+1)=max(audioPwr(:,k)*alpha + activity(:,k)*(1-alpha),activity(:,k));
    end
    audioPwr(:,1)=audioPwr(:,end);    
    
    % averaged energy
    energy=sum(audioPwr,2)/mAvg;
    % spectrum of audio from average window
    spect=(MneurToBin*energy).*exp(1i*phs);
    
    % reduce spectrum to tone frequency components
    toneMags = interp1(fftFreqs,abs(spect),toneFreqs,'linear','extrap');
    tonePhases = interp1(fftFreqs,angle(spect),toneFreqs,'linear','extrap');
    
    interpSpect(:,blkNumber) = toneMags.*exp(1j*tonePhases);    
end

    %interpolated spectral envelope scaling
    specVec = (1:blkNumber)*nFFT/2;
    newTimeVec = 0:(nBlocks-(nFFT/2));
    interpSpect2 = zeros(length(toneFreqs),length(newTimeVec));
    modTones = zeros(size(interpSpect2));
    
    %interpolate spectral envelope from frames back to td
    for freq = 1: length(toneFreqs)
        tEnvMag = interp1(specVec,abs(interpSpect(freq,:)),newTimeVec,'linear','extrap');
        tEnvPhs = interp1(specVec,angle(interpSpect(freq,:)),newTimeVec,'linear','extrap');
        interpSpect2(freq,:) = tEnvMag.*exp(1j*tEnvPhs);
        modTones(freq,:) = tones(freq,1:end-(nFFT/2-1)).*abs(interpSpect2(freq,:));
    end
    
    audioOut = sum(modTones,1);

end

