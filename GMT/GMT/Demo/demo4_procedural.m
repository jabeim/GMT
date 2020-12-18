% This demo implements the exact same processing strategy as
% demo3_SpecRes_NR.m in a purely procedural style, i.e. avoiding
% use of the object-oriented framework built around the core functionality
% of each processing block. Corresponding lines from the object
% oriented demo are quoted in comments.

%% Global strategy parameters 
par_strat = struct( ...
    'fs', 17400, ...
    'nFft', 256, ...
    'nHop', 20,  ...
    'nChan', 15,  ...
    'startBin', 6,  ...
    'nBinLims', [2, 2, 1, 2, 2, 2, 3, 4, 4, 5, 6, 7, 8, 10, 56], ...
    'window', 0.5*(blackman(256) + hanning(256)), ...
    'pulseWidth', 18, ...
    'verbose', 0, ...
    'wavFile', ['Sounds' filesep 'AzBio_3sent_65dBSPL.wav']...
    );

%% Setting up parameter structures for each subsequent function call

% src = ReadWavUnit(strat, 'SRC', 'Sounds\AzBio_3sent.wav');
par_readWav = struct( ...
    'parent', par_strat, ...
    'tStartEnd', [], ... 
    'iChannel', 1 ...
    );

% mix = AudioMixerUnit(strat, 'MIX', 1, 0, 'rel');   % 1 input, no re-scaling, assuming correct scaling in input wav file (parameters: 0 dB rel. to input)
% (.. single input, transparent gain -> omitted ..)

% pre = HarmonyPreemphasisUnit(strat, 'PRE');             % pre-emphasis filter
par_pre = struct(... 
   'parent', par_strat, ...
   'coeffNum', [0.7688 -1.5376 0.7688], ...  % numerator coefficients
   'coeffDenom', [1 -1.5299 0.5453]     ...  % denominator coefficients
   );

% agc = DualLoopTdAgcUnit(strat, 'AGC');                  % AGC
par_agc = struct(...
    'parent', par_strat, ...
    'kneePt', 4.476, ... % compression threshold [log2]
    'compRatio', 12, ... % compression ratio above knee-point (in log-log space) [> 1] [12]
    'tauRelFast', -8 / (17400 * log(0.9901)) * 1000, ... % fast release time const [ms] [46.21]
    'tauAttFast', -8 / (17400 * log(0.25)) * 1000, ...   % fast attack time const  [ms] [0.33]
    'tauRelSlow', -8 / (17400 * log(0.9988)) * 1000, ... % slow release time const [ms] [382.91]
    'tauAttSlow', -8 / (17400 * log(0.9967)) * 1000, ... % slow attack time const  [ms] [139.09]
    'maxHold', 1305, ...    % max. hold counter value [int >= 0] [1305]
    'g0', 6.908, ...        % gain for levels < kneepoint [log2] [6.908; approx = 41.6dB]
    'fastThreshRel', 8, ... % relative threshold for fast loop [dB] [8]
    'cSlowInit', 0.5e-3, ...   % initial value for slow averager, 0..1, [] for auto; default: 1 
    'cFastInit', 0.5e-3, ...   % initial value for fast averager, 0..1, [] for auto; default: 1 
    'controlMode', 'naida', ... % how to use control signal, if provided on port #2? ['naida' / 'direct'] ['naida']
    'clipMode', 'limit', ... % output clipping behavio ['none' / 'limit' / 'overflow'] ['none']
    'decFact', 8, ...       % decimation factor (i.e. frame advance)
    'envBufLen', 32, ...    % buffer (i.e. frame) length for envelope computation
    'gainBufLen', 16, ...   % buffer length for gain smoothing
    'envCoefs', [-19,  55,   153,  277,  426,  596,	784,  983,   ...  % tapered envelope data window  
                 1189, 1393, 1587, 1763,	1915, 2035,	2118, 2160,  ...
                 2160, 2118,	2035, 1915, 1763, 1587,	1393, 1189,  ...
                 983,  784,	596,  426,	277,  153,  55,   -19 ] / (2^16) ...
    );

% wb = WinBufUnit(strat, 'WB');                           % buffering and windowing
par_winBuf = struct( ...
    'parent', par_strat, ...
    'bufOpt', []  ... 
    );

% fftfb = FftFilterbankUnit(strat, 'FFT');                % FFT 
par_fft = struct( ...
    'parent', par_strat, ...
    'combineDcNy', false, ...           % Combine DC and Nyquist bins into single complex 1st bin?
    'compensateFftLength', false, ...   % Divide FFT coefficients by nFft/2? [boolean]
    'includeNyquistBin', false ...      % Return bin #nFft/2+1 in output? [boolean] 
    );

% env = HilbertEnvelopeUnit(strat, 'HILB');               % Hilbert envelopes 
par_hilbert = struct( ...
    'parent', par_strat, ...
    'outputOffset', 0, ...      % scalar offset added to all channel outputs; [log2] [0] Use with caution!
    'outputLowerBound', 0, ...  % lower bound applied to output (after offset) [log2] [0]
    'outputUpperBound', Inf ... % lower bound applied to output (after offset) [log2] [Inf]
    );

% engy = ChannelEnergyUnit(strat, 'ENGY', 2);             % channel energies (for noise reduction SNR estimation); 2 inputs (to account for AGC gain)
par_energy = struct( ...
    'parent', par_strat, ...
    'gainDomain', 'linear' ...  % domain of gain input (#2)  ['linear','db','log2']
    );

% nr = NoiseReductionUnit(strat, 'NR', 1, 'log2', false);     % noise reduction; 'log2' makes gain output commensurable with Hilbert envelopes
par_nr = struct( ...
    'parent', par_strat,  ...
    'gainDomain', 'log2', ...   % domain of gain output on port 2 (if applicable) ['linear','db','log2'] ['linear']
    'tau_speech', 0.0258, ...   % time constant of speech estimator [s] [0.0258]
    'tau_noise', 0.219, ...     % time constant of noise estimator [s] [0.219]
    'threshHold', 3, ...        % hold threshold (onset detection criterion) [dB, > 0] [3]
    'durHold', 1.6, ...         % hold duration (following onset) [s] [1.6]
    'maxAtt',  -12, ...         % maximum attenuation (applied for SNRs <= snrFloor) [dB] [-12]
    'snrFloor', -2, ...         % SNR below which the attenuation is clipped [dB] [-2]
    'snrCeil', 45, ...          % SNR above which the gain is clipped  [dB] [45]
    'snrSlope', 6.5, ...        % SNR at which gain curve is steepest  [dB] [6.5]
    'slopeFact', 0.2, ...       % factor determining the steepness of the gain curve [> 0] [0.2]
    'noiseEstDecimation', 1, ...    % down-sampling factor (re. frame rate) for noise estimate [int > 0] [1]
    'enableContinuous', false, ...  % save/restore states across repeated calls of run [bool] [false]
    'initState', struct('V_s', -30 * ones(15,1), 'V_n', -30 * ones(15,1)) ...    % initial state
    );

% gapp = ElementwiseUnit(strat, 'GAPP', 2, @plus, true);  % NR gain application: element-by-element sum of 2 input;
%  (no corresponding parameter struct, addition handled in the function call section below)

% spl = SpecPeakLocatorUnit(strat, 'SPL');                % channel peak frequency and target location estimation
par_peak = struct( ...
    'parent', par_strat, ...
    'binToLocMap', [zeros(1,6), 256, 640, 896, 1280, 1664, 1920, 2176, ...   % 1 x nBin vector of nominal cochlear locations for the center frequencies of each STFT bin
                    2432, 2688, 2944, 3157, 3328, 3499, 3648, 3776, 3904, 4032, ...   % as in firmware; values from 0 .. 15 (originally in Q9 format)
                    4160, 4288, 4416, 4544, 4659, 4762, 4864, 4966, 5069, 5163, ...   % corresponding to the nominal steering location for each 
                    5248, 5333, 5419, 5504, 5589, 5669, 5742, 5815, 5888, 5961, ...   % FFT bin
                    6034, 6107, 6176, 6240, 6304, 6368, 6432, 6496, 6560, 6624, ...
                    6682, 6733, 6784, 6835, 6886, 6938, 6989, 7040, 7091, 7142, ...
                    7189, 7232, 7275, 7317, 7360, 7403, 7445, 7488, 7531, 7573, ...
                    7616, 7659, 7679 * ones(1,53)] / 512 ...
    );

% csw = CurrentSteeringWeightsUnit(strat, 'CSW');         % current steering weights based on target location
par_steer = struct( ...
    'parent', par_strat, ...
    'nDiscreteSteps', 9, ...   % nr. of discretization steps  [int >= 0] [9]; 0 -> no discretization
    'steeringRange', 1.0 ...  % steering range between electrodes [0..1] [1.0]    
    );

% csynth = CarrierSynthesisUnit(strat, 'CSYNTH');         % synthesize electrode carrier signal at FT rate (temporal fine structure)
par_carrierSynth = struct( ...
    'parent', par_strat, ...
    'fModOn', 0.5,  ... % peak frequency up to which max. modulation depth is applied [fraction of FT rate] [0.5]
    'fModOff', 1.0, ... % peak frequency beyond which no modulation is applied  [fraction of FT rate] [1.0]
    'maxModDepth', 1.0, ... % maximum modulation depth [0.0 .. 1.0] [1.0]
    'deltaPhaseMax', 0.5 ...% maximum phase rotation per FT frame [turns] [0.5] 
    );                      % Set <= 0.5 to avoid aliasing for fPeak > FT_rate/2    

% map = F120MappingUnit(strat, 'MAP');                    % combine envelopes and carriers, and map to stimulation current amplitude 
par_mapper = struct( ...
    'parent', par_strat, ...
    'mapM', 500 * ones(1,16), ...       % M levels [uAmp] 
    'mapT', 50  * ones(1,16), ...       % T levels [uAmp] 
    'mapIdr', 60 * ones(1,16), ...      % IDRs  [dB]
    'mapGain', 0 * ones(1,16), ...      % channel gains [dB]
    'mapClip', 2048 * ones(1,16), ...   % clipping level [uAmp] [2048]
    'chanToElecPair', 1:15, ...         % 1 x nChan vector defining mapping of logical channels to electrode pairs (1 = E1/E2, ...) [] [1:nChan]
    'carrierMode', 1 ...                % carrierMode - how to apply carrier [0 - no carrier, 1 - to input, 2 - to output] [1]
    );

% egram = F120ElectrodogramUnit(strat, 'EGRAM', true);
par_elgram = struct(...
    'parent', par_strat, ...
    'cathodicFirst', true, ...  %  start pulse with cathodic phase? [bool] 
    'channelOrder', [1 5 9 13 2 6 10 14 3 7 11 15 4 8 12], ... % default F120 staggering order [DO NOT MODIFY]
    'enablePlot', true, ...     %  plot electrodogram? [bool]
    'colorScheme',2,...
    'outputFs', [] ...     %  output sampling frequency [Hz]. LEAVE EMPTY TO DEFAULT TO NATIVE RATE (55555.55...Hz)
    );

par_validate = struct(...
    'parent',par_strat,...
    'lengthTolerance',50,...                % maximum allowable difference between validation file and submitted data
    'saveIfSimilar',true,...                % force saving files even if result is similar to default processing [FOR TESTING]
    'differenceThreshold',1,...             % minimum current difference across time for each channel
    'maxSimilarChannels',8,...              % max number of similar channels allowed before preventing save
    'elgramFs',par_elgram.outputFs,...      % match electrodogram sampling frequency
    'outFile',''...
    );

par_vocoder = struct(...
    'parent',par_strat,...,
    'audioFs',48000,...
    'saveAudioOutput',true,...
    'audioOutputFile','' ...
    );
%% Function calls
% ( strat.connect(block_1, block_2), ..., strat.run() ) 

sig_smp_wavIn                           = readWavFunc(par_readWav); % read wav input; assumes correct scaling of 111.6 dB SPL peak full-scale in wav input

sig_smp_wavPre                          = tdFilterFunc(par_pre, sig_smp_wavIn); % pre-emphasis
[sig_smp_wavAgc, sig_smp_gainAgc]       = dualLoopTdAgcFunc(par_agc, sig_smp_wavPre); % AGC
sig_frm_audBuffers                      = winBufFunc(par_winBuf, sig_smp_wavAgc); % buffering 
sig_frm_fft                             = fftFilterbankFunc(par_fft, sig_frm_audBuffers); % STFT

sig_frm_hilbert                         = hilbertEnvelopeFunc(par_hilbert, sig_frm_fft); % Hilbert envelopes
sig_frm_energy                          = channelEnergyFunc(par_energy, sig_frm_fft, sig_smp_gainAgc); % channel energy estimates
sig_frm_gainNr                          = noiseReductionFunc(par_nr, sig_frm_energy); % noise reduction
sig_frm_hilbertMod                      = sig_frm_hilbert + sig_frm_gainNr; % apply noise reduction gains to envelopes

%   sub-sample every third FFT input frame
sig_3frm_fft = sig_frm_fft(:,3:3:end);    
[sig_3frm_peakFreq, sig_3frm_peakLoc]   = specPeakLocatorFunc(par_peak, sig_3frm_fft); % peak frequency and location estimation 
%   up-sample back to full frame rate (and add padding for skipped initial frames)
sig_frm_peakFreq = repelem(sig_3frm_peakFreq,1,3);                         
sig_frm_peakFreq = [zeros(size(sig_frm_peakFreq, 1),2), sig_frm_peakFreq]; 
sig_frm_peakFreq = sig_frm_peakFreq(:,1:size(sig_frm_fft, 2)); 
sig_frm_peakLoc = repelem(sig_3frm_peakLoc,1,3); 
sig_frm_peakLoc = [zeros(size(sig_frm_peakLoc, 1),2), sig_frm_peakLoc];
sig_frm_peakLoc = sig_frm_peakLoc(:,1:size(sig_frm_fft, 2)); 

sig_frm_steerWeights                    = currentSteeringWeightsFunc(par_steer, sig_frm_peakLoc);  % current steering, based on peak location
[sig_ft_carrier, sig_ft_idxFtToFrm]     = carrierSynthesisFunc(par_carrierSynth, sig_frm_peakFreq); % carrier synthesis, based on peak frequencies
sig_ft_ampWords                         = f120MappingFunc(par_mapper, sig_ft_carrier, ... % combine envelepes, carrier and current steering weights, compute current outputs
                                            sig_frm_hilbertMod, sig_frm_steerWeights, sig_ft_idxFtToFrm); 

elgram = f120ElectrodogramFunc(par_elgram, sig_ft_ampWords);  % plot electrogram
saved = validateOutputFunc(par_validate,elgram);

[audioOut,audioFs] = vocoderFunc(par_vocoder,elgram);

sound(audioOut,audioFs);
