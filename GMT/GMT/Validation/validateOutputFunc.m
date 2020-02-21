function saved = validateOutputFunc(par,electrodogram)
    if isstr(electrodogram)
        loaded = load(electrodogram);
        electrodogram = loaded.elData;
    end
    
    if par.elgramFs ~= 200e3; error('Electrodogram must be generated with 200 kHz sampling rate'); end
        
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
  
    validationFileName = [inputFileName(fnStart(end)+1:strfind(inputFileName,'.wav')-1) '_validation.mat'];
    
    defaultData = load(['Validation\' validationFileName]);
    outputDifference = sum(electrodogram-full(defaultData.elData),2);
    
    if par.saveWithoutValidation == true
        if any(outputDifference < par.differenceThreshold)
            channels = find(outputDifference < par.differenceThreshold)';
            if length(channels) == 1
                disp(['Channel: ' num2str(channels) ' is too similar to the default output.'])
            else
                disp(['Channels: ' num2str(channels) ', are too similar to the default output.'])
            end
        end
        elData = sparse([outputDifference electrodogram]);
        
        if length(par.outFile) == 0
            timestr = datestr(now,'yyyymmdd_HHMMSS');
            save(['Output/elGramOutput_' timestr '.mat'], 'elData')
        else
            save(['Output/' par.outFile],'elData')
        end
        saved = true;
        
    else
        if any(outputDifference < par.differenceThreshold)
            channels = find(outputDifference < par.differenceThreshold)';
            if length(channels) > par.maxSimilarChannels;
                if length(channels) == 1
                    disp(['Channel: ' num2str(channels) ' is too similar to the default output. DATA NOT SAVED!'])
                else
                    disp(['Channels: ' num2str(channels) ', are too similar to the default output. DATA NOT SAVED'])
                end
                saved = false;
            end
        else
            elData = sparse([outputDifference electrodogram]);
            if length(par.outFile) == 0
                timestr = datestr(now,'yyyymmdd_HHMMSS');
                save(['Output/elGramOutput_' timestr '.mat'], 'elData')
            else
                save(['Output/' par.outFile],'elData')
            end
            saved = true;
        end     
    end
    