function [audioOut,audioFs] = vocoderFunc(par,electrodogram)
% scan input parameters adding necessary default values
defaultParNames = {'captFs','nCarriers','elecFreqs','spread','neuralLocsOct','nNeuralLocs','MCLmuA','TmuA','tAvg','audioFs','tPlay','tauEnvMS','nl','resistorVal'};
defaultParValues = {200000,20,[],[],[],100,[],[],.005,48000,[],10,5,2};

for i = 1:length(defaultParNames)
    if ~isfield(par,defaultParNames{i})
        par.(defaultParNames{i}) = defaultParValues{i};
    end
end

%% scale and preprocess electrodogram data
scale2MuA = 1000/par.resistorVal;
electrodeAmp = electrodogram;
nElec = size(electrodogram,1);
elData = electrodeAmp*scale2MuA;
captTs = 1/captFs;


%% compute electrode locations in terms of frequency
if isesmpty(par.elecFreqs)
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
alpha = exp(-1/(taus*captFs));

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
    TmuA = 500*ones(1,nElec);
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

%% length of audio
if isempty(par.tPlay)
    tPlay = 10*tWin;
else
    tPlay = ceil(par.tPlay/tWin)*twin;
end

%charge 2 EF matrix
charge2EF = zeros(nNeuralLocs,nElec);
elecFreqOct = log2(elecFreqs);

for iEl = 1:nElec  
    x = spread(elecPlacement(iEl)).fOct+elecFreqOct(iEl);
    xx = spread(elecPlacement(iEl)).amp;
    xxx = neuralLocsOct;
    
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
%MneurToBin(:,end)=0;
%imagesc(MneurToBin);
%win=Leos_window(nFFT)';
%% window shape |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
win=.5-.5*cos(2*pi*(0:(nFFT-1))/(nFFT-1));
%win=.5*ones(size(win));

%% other auxiliary variables
playOverAvgRatio=round(tPlay/tAvg);
stateHolder=zeros(nFFT+(playOverAvgRatio-1)*nAvg,1);
shli=1;%state holder lower index


phs=2*pi*rand(floor(nFFT/2),1);%|||||||||||||||||||||||||||||||||||||||||||
%phs=zeros(floor(nFFT/2),1);
%dphi=(pi*nAvg/floor(nFFT/2))*(1:floor(nFFT/2))';
dphi=2*pi*(1:floor(nFFT/2))'*nAvg/nFFT; %||||||||||||||||||||||||||||||||||

audioPwr=zeros(nNeuralLocs,blkSize+1);

%actual iteration over blocks of current
% normRamp=charge2EF*unipolarPulseWidthMuSec/(MCLnC-TnC);
% normOffset=(TnC)/(MCLnC-TnC);
%whos charge2EF
%whos MCLnC
 
%normRamp=bsxfun(@times,charge2EF*unipolarPulseWidthMuSec,1./(MCLmuA(:)'-TmuA));
%normOffset=charge2EF*((TmuA(:))./(MCLmuA(:)-TmuA(:)));
M=interp1(elecFreqOct,MCLmuA,neuralLocsOct);M=M(:);
M(neuralLocsOct<elecFreqOct(1))=MCLmuA(1);
M(neuralLocsOct>elecFreqOct(nElec))=MCLmuA(nElec);
T=interp1(elecFreqOct,TmuA,neuralLocsOct);T=T(:);
T(neuralLocsOct<elecFreqOct(1))=TmuA(1);
T(neuralLocsOct>elecFreqOct(nElec))=TmuA(nElec);
normRamp=bsxfun(@times,charge2EF,1./(M-T));
normOffset=T./(M-T);

elData(elData<0)=0;%if abs is used, ireelevant electrodes add constructively
fftFreqs=(1:floor(nFFT/2))*audioFs/nFFT;
fftFreqsOct=log2(fftFreqs);


% generate tone complex
nBlocks = (nFFT/2*floor(size(eldata,2)/blkSize+1))-1;

tones = zeros(nBlocks,nCarriers);
toneFreqs = generate_cfs(20,20000,nCarriers);
t =(1:nBlocks)/audioFs;

for toneNum = 1:nCarriers
    tones(:,toneNum) = sin(2*pi*toneFreqs(toneNum)*t+phs(toneNum));
end

interpSpect = zeros(nCarriers,floor(size(elData,2)/blkSize));
fftFreqs = 1:floor(nFFT/2)+1*audioFs/nFFT;


%audioFs/nFFT
%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
for blkNumber=(1):floor(size(elData,2)/blkSize)
    %% charge to patient-normalized electric field
    timeIdx=((blkNumber-1)*blkSize+1) : (blkNumber*blkSize);
    
   
    %EF=normRamp*abs(elData(:,timeIdx))-normOffset;
    EF=max(0,bsxfun(@plus,normRamp*elData(:,timeIdx),-normOffset));
    
    
    % normalized EF to activity
    nl = par.nl;
    activity=max(0,min(exp(-nl+nl*EF),1)-exp(-nl))/(1-exp(-nl));%||||||||||||||||||||||||||||||||||||||||  
    
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
    
    interpSpect(:,blkNumber) = toneMags*exp(1j*tonePhases);
    
end

    %interpolated spectral envelope scaling
    
    specVec = 1:blkNumber*nFFT/2;
    newTimeVec = 1:(nBlocks-(nFFT/2));
    interpSpect2 = zeros(length(toneFreqs),length(newTImeVec));
    modTones = zeros(size(interpSpect2));
    
    %interpolate spectral envelope from frames back to td
    for freq = 1: length(toneFreqs)
        tEnvMag = interp1(specVec,abs(interpSpect(freq,:)),newTimeVec,'linear','extrap');
        tEnvPhs = interp1(specVec,angle(interpSpect(freq,:)),newTimeVec,'linear','extrap');
        interpSpect2(freq,:) = tEnvMag*exp(1j*tEnvPhs);
        modTones(freq,:) = tones(1:end-nFFT/2,freq)*abs(interpSpect2(freq,:));
    end
    
    audioOut = sum(modTones,1);

end

