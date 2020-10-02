function saved = validateOutputFunc(par,electrodogram)
% saved = validateOutputFunc(par,electrodogram)
%
% Validates electrodogram outputs based on contest rules by ensuring that:
% 1. electrodogram sampling rate is 55556 Hz
% 2. electrodogram contains 16 rows (channels)
% 3. electrodogram length matches the resampled source audio (within a small tolerance)
% 4. electrodogram channels are charge-balanced: abs(sum(channel)) < epsilon
%
% Validation also attempts to load a validation file containing the
% electrodogram of the same stimulus preprocessed by the default strategy
% for comparison. A warning is issued if no validation file is found, and
% the similarity comparison process is skipped.
%
% A simple subtractive analysis of each channel is used to estimate the
% similarity between elecdtrode channels using cross correllation to time
% align the data.
%
% INPUT:
%   electrodogram - either a 16 x n matrix, where N is audioDuration*55556
%           samples, or a path to a .mat file containing the electrodogram data
%           saved as a variable 'elData'
%
% FIELDS FOR PAR:
%   parent.wavFile - the name of the source audio file being processed, ['Sounds\example.wav']
%   lengthTolerance - the number of samples difference allowable between validation and electrodogram array lengths, [15]
%   saveIfSimilar - whether or not to save the output matrix if similarity between electrodogram and default data is high, [bool]
%   differenceThreshold - max value for sum(abs(electrodogram-validationData)) in each channel, 
%                         channels are flagged as similar if the result is less than the threshold, [int]
%   maxSimilarChannels - maximum number of channels where the similarity exceeds the difference Threshold, [8]
%   elGramFs - the sampling frequency used to generate electrodogram; also stored in the parameters for f120ElectrodogramFunc
%               MUST BE SET TO 55556 Hz
%   outFile - additional text appended to the filename of electrodogram data saved by the validation function, '' [string]
%
% OUTPUT:
%    saved - boolean indicating whether or not the data in electrodogram was saved to an output file, true [bool]
%
    % assume that
    inputDataIsMatrix = true;

    if isstr(electrodogram)
        loaded = load(electrodogram);
        electrodogram = loaded.elData;
        inputDataIsMatrix = false;
    end
    
    if par.elgramFs ~= 55556; error('Electrodogram must be generated with 55556 Hz sampling rate'); end        
    if length(size(electrodogram)) ~= 2; error('Electrodogram must be a 2 dimensional matrix!'); end
    if size(electrodogram,1) ~= 16
        if size(electrodogram,2) ~=16
            error(['Electrodogram dimensions should be 16 x numSamples, currently ' num2str(size(electrodogram,1)) ' x ' num2str(size(electrodogram,2)) '.'])
        else
            electrodogram = electrodogram';
        end
    end
    
    inputFileName = par.parent.wavFile;
    fnStart = strfind(inputFileName,filesep);
    validationFileName = inputFileName(fnStart(end)+1:strfind(inputFileName,'.wav')-1);
    skipMatrixSubtraction = false;
    
    try
        defaultData = load(['Validation' filesep validationFileName '_validation.mat']);
        validationData = full(defaultData.elData); % validation files are sparse matrices by default, conversion to full is necessary for xcorr
    catch err
        % add some verbosity about file not found
        if err.identifier == 'MATLAB:load:couldNotReadFile'
            warning(['No validation file found for source file: ' inputFileName ' skipping matrix similarity comparison!'])
            disp('Ensure validation files are stored in the validation folder. Visit our website to download the validation files if you do not have them.')
            audioSourceInfo = audioinfo(inputFileName);
            validationData = zeros(16,fix(audioSourceInfo.TotalSamples/audioSourceInfo.SampleRate*55556));  % just a matrix of zeros for size comparisons
            skipMatrixSubtraction =true;
        else
            rethrow(err)
        end
    end
   
    % compare length to validation file
    if skipMatrixSubtraction
        % if there was no validation file, skip this step and continue
        outputDifference = [];  % leave empty validation matrix, this will pass through validation as a case where no channels are similar
    else
        if size(validationData,2)-par.lengthTolerance <= size(electrodogram,2) < size(validationData,2)+par.lengthTolerance
            for i = 1:size(validationData,1)
                outputDifference(1,i) = xCorrSimilarity(validationData(i,:),full(electrodogram(i,:)));
            end
        else
            error(['Electrodogram length exceeds validation tolerance. Expected: ' num2str(size(validationData,2)) ' samples. +-' num2str(par.lengthTolerance) '%, found: ' num2str(size(electrodogram,2)) 'samples.'])
        end
    end
    
    
    % check for charge balancing
    chargeBalance = zeros(size(electrodogram,1));
    for i = 1:size(electrodogram,1)
        if sum(electrodogram(i,:)) > eps
            chargeBalance(i) = 1;
        end
    end
   
    if sum(chargeBalance) > 0
        error(['Channels: ' numstsr(find(chargeBalance)) ' are not charge balanced. Within-channel current must sum to 0.'])
    end
    

    elData = sparse(electrodogram);
    
    % file saving logic depending on whether or not data or an existing
    % file was passed to the function
    if inputDataIsMatrix == true
        % allow saving output matrix regardless of similarity, but still
        % print similarity info to console
        if par.saveIfSimilar == true
            if any(outputDifference < par.differenceThreshold)
                channels = find(outputDifference < par.differenceThreshold);
                if length(channels) == 1
                    disp(['Channel: ' num2str(channels) ' is very similar to the default output.'])
                elseif length(channels) > 1
                    disp(['Channels: ' num2str(channels) ', are very similar to the default output.'])
                end
            end
            
            if length(par.outFile) == 0
                % generate a timestamp if user provides no other description in par.outFile
                timestr = datestr(now,'yyyymmdd_HHMMSS');
                save(['Output' filesep validationFileName '_elGramOutput_' timestr '.mat'], 'elData')
            else
                save(['Output' filesep validationFileName '_elGramOutput_' par.outFile],'elData')
            end
            saved = true;
        % if similarity is used to prevent saving, follow this logic. Use
        % warnings to call attention to lack of data saving.
         else
            channels = find(outputDifference < par.differenceThreshold);
            if length(channels) > par.maxSimilarChannels
                % display warnings if similarity validation fails
                if length(channels) == 1
                    warning(['Channel: ' num2str(channels) ' is very similar to the default output. DATA NOT SAVED!'])
                elseif length(channels) > 1
                    warning(['Channels: ' num2str(channels) ', are very similar to the default output. DATA NOT SAVED'])
                end
                saved = false;
            else
                % save files if similarity validation passes.
                if length(par.outFile) == 0
                    % generate a timestamp if user provides no other description in par.outFile
                    timestr = datestr(now,'yyyymmdd_HHMMSS');
                    save(['Output' filesep validationFileName '_elGramOutput_' timestr '.mat'], 'elData')
                else
                    save(['Output' filesep validationFileName '_elGramOutput_' par.outFile],'elData')
                end
                saved = true;
            end      
        end
    end
    