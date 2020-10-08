% This demo implements SpecRes with noise reduction and computes the electrode
% output for some speech tokens, as measured using a load-board and scope
initGmtClassPath;

%% Create strategy
strat = FftStrategy(); % AzBio_3sent.wav is assigned as default sound input, use line below to load other files
strat.wavFile = 'AzBio_3sent_65dBSPL.wav';  % same as default, but scaled to 65dB SPL RMS assuming 111.6dB 
%% Create instances of ProcUnits and add them to strategy
src = ReadWavUnit(strat, 'SRC');  % use wav file as input
mix = AudioMixerUnit(strat, 'MIX', 1, 0, 'rel');        % 1 input, no re-scaling, assuming correct scaling in input wav file (parameters: 0 dB rel. to input)

pre = HarmonyPreemphasisUnit(strat, 'PRE');             % pre-emphasis filter
agc = DualLoopTdAgcUnit(strat, 'AGC');                  % AGC 
wb = WinBufUnit(strat, 'WB');                           % buffering and windowing
fftfb = FftFilterbankUnit(strat, 'FFT');                % FFT 

env = HilbertEnvelopeUnit(strat, 'HILB');               % Hilbert envelopes 
engy = ChannelEnergyUnit(strat, 'ENGY', 2);             % channel energies (for noise reduction SNR estimation); 2 inputs (to account for AGC gain)
nr = NoiseReductionUnit(strat, 'NR', 1, 'log2', false); % noise reduction; 'log2' makes gain output commensurable with Hilbert envelopes
gapp = ElementwiseUnit(strat, 'GAPP', 2, @plus, true);  % NR gain application: element-by-element sum of 2 input;
                                                        % ('true' indicates that @add supports matrix inputs natively)
                                                        
spl = SpecPeakLocatorUnit(strat, 'SPL');                % channel peak frequency and target location estimation
csw = CurrentSteeringWeightsUnit(strat, 'CSW');         % current steering weights based on target location
csynth = CarrierSynthesisUnit(strat, 'CSYNTH');         % synthesize electrode carrier signal at FT rate (temporal fine structure)

map = F120MappingUnit(strat, 'MAP');                    % combine envelopes and carriers, and map to stimulation current amplitude 

egram = F120ElectrodogramUnit(strat, 'EGRAM', true);    % generate electrodogram from mapper output, plotting enabled


%% set (non-default) block parameters
agc.cFastInit = 0.5e-3; 
agc.cSlowInit = 0.5e-3; % start AGC in fully "relaxed" state
agc.clipMode = 'limit';

nr.initState = struct('V_s', -30 * ones(15,1), 'V_n', -30 * ones(15,1));

map.mapM = 500 * ones(1,16); % M = 500 uA
map.mapT = 50  * ones(1,16); % T = 50 uA


%% connect ProcUnits (using block labels)
strat.connect(src, mix);
strat.connect(mix, pre);

strat.connect(pre, agc);
strat.connect(agc, wb);
strat.connect(wb, fftfb);  

strat.connect(fftfb, env);
strat.connect(fftfb, engy);
strat.connect(agc, 2, engy, 2);
strat.connect(engy, nr);
strat.connect(env, gapp);
strat.connect(nr, gapp, 2);

strat.connect(fftfb, spl);
strat.connect(spl, csynth);
strat.connect(spl, 2, csw);

strat.connect(gapp, map); 
strat.connect(csw, map, 2);
strat.connect(csynth, map, 3);
strat.connect(csynth, 2, map, 4);

strat.connect(map, egram);

egram.outputFs = 55556;   % 200 kHz scope sampling rate 
egram.resistance = 10000;  % 10 kOhm load-board resistors
egram.colorScheme = 4;

% create "offline" viewer (display strategy as is, no dynamic updating)
hFig = csViewer(strat, [], 0); 

%% run strategy
strat.run();
