function saved = validateOutputFunc(par,electrodogram)
    if isstr(electrodogram)
        loaded = load(electrodogram);
        electrodogram = loaded.elData;
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
    
    try
        defaultData = load(['Validation' filesep validationFileName '_validation.mat']);
    catch err
        % add some verbosity about file not found
        if err
            if par.skipValidation == true
                warning(['No validation file found for source file: ' inputFileName ' skipping matrix similarity comparison!'])
                defaultData.elData = electrodogram;
            else
                rethrow(err);
            end

        end
    end
   
    % compare length to validation file
    if size(defaultData.elData,2)*.99 <= size(electrodogram,2) < size(defaultData.elData,2)*1.01
        
    else
        error(['Electrodogram length exceeds validation tolerance. Expected: ' num2str(size(defaultData,2)) ' samples. +-' num2str(lengthTolerance) '%, found: ' num2str(size(electrodogram,2)) 'samples.'])
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
    
    % matrix subtraction comparison to catch exact duplicates
    if size(electrodogram) == size(defaultData,elData)
        outputDifference = sum(electrodogram-full(defaultData.elData),2);
    else
        outputDifference = zeros(size(electrodogram,1));
    end
    elData = sparse(electrodogram);
    
    if par.saveWithoutValidation == true
        if any(outputDifference < par.differenceThreshold)
            channels = find(outputDifference < par.differenceThreshold)';
            if length(channels) == 1
                disp(['Channel: ' num2str(channels) ' is too similar to the default output.'])
            else
                disp(['Channels: ' num2str(channels) ', are too similar to the default output.'])
            end
        end
%         elData = sparse([outputDifference electrodogram]);
        
        if length(par.outFile) == 0
            timestr = datestr(now,'yyyymmdd_HHMMSS');
            save(['Output/elGramOutput_' timestr '.mat'], 'elData')
%             csvwrite(['Output/elGramOutput_' timestr '.dat'], elData)
        else
            save(['Output/' par.outFile],'elData')
%             csvwrite(['Output/' par.outFile],elData)
        end
        saved = true;
        
    else
        channels = find(outputDifference < par.differenceThreshold)';
        if length(channels) > par.maxSimilarChannels;
            if length(channels) == 1
                disp(['Channel: ' num2str(channels) ' is too similar to the default output. DATA NOT SAVED!'])
            else
                disp(['Channels: ' num2str(channels) ', are too similar to the default output. DATA NOT SAVED'])
            end
            saved = false;
        else
%             elData = sparse([outputDifference electrodogram]);
            if length(par.outFile) == 0
                timestr = datestr(now,'yyyymmdd_HHMMSS');
                save(['Output' filesep validationFileName '_elGramOutput_' timestr '.mat'], 'elData')
%                 csvwrite(['Output/elGramOutput_' timestr '.dat'], elData)
            else
                save(['Output' filesep validationFileName '_elGramOutput_' par.outFile],'elData')
%                 csvwrite(['Output/' par.outFile],elData)
            end
            saved = true;
        end      
    end
    