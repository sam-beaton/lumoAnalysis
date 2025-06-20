
function [parcelBlockAveraged, params] = getParcelAverageBlock(parcelAveraged, info, params)

% Calculates the parcel-wise block average of the task data for each
% chromophore, using the parcel-wise average (parcelAveraged) full time
% series data

    parcelBlockAveraged = cell(info.io.Nwl,1); %should be 2*1
    
    % account for different sampling frequencies
    dtPre = floor(params.dtPre*(info.system.framerate/10)); 
    dtAfter = floor(params.dtAfter*(info.system.framerate/10));

    % get pulse indices
    allPulses = [];
    % loop through Pulse_2 to Pulse_6 (Pulse_1 is baseline)
    for i = 2:6
        fieldName = sprintf('Pulse_%d', i);
        
        % if the field exists
        if isfield(info.paradigmFull, fieldName)
            thisPulse = info.paradigmFull.(fieldName); %give temp name          
            allPulses = [allPulses; thisPulse(:)]; % add to existing
        end
    end
    % get sample point where pulse/stim marker occurs
    allPulses = info.paradigmFull.synchpts(allPulses);
    [allPulses, params.blockRemoved] = analysisTools.chopBlocks(allPulses, size(parcelAveraged{2}, 2), dtPre, dtAfter);
    
    % store block data
    % get number of pulses
    numBlocks = length(allPulses);
    numSamples = dtPre + dtAfter + 1;
    timeAxis = linspace(-params.dtPre, params.dtAfter, numSamples);
    
    for iChrom = 1:info.io.Nwl
    
        parcelData = parcelAveraged{iChrom};
        
        %initialise storage for this chromophore
        numParcels = size(parcelAveraged{iChrom}, 1);
        parcelBlockAveragedChrom = zeros(numParcels, numBlocks, (dtAfter+1));
    
        for iParc = 1:numParcels
            for iBlock = 1:numBlocks
                %unnecessary after removing pulses in chopBlocks but kept
                %in for second check layer:
                if allPulses(iBlock)-dtPre >= 1 && allPulses(iBlock)+dtAfter <= size(parcelData, 2) 
                    parcelBlockAveragedChrom(iParc, iBlock, :) = parcelData(iParc, allPulses(iBlock):allPulses(iBlock)+dtAfter);
                    parcelBlockAveragedChrom(iParc, iBlock, :) = parcelBlockAveragedChrom(iParc, iBlock, :) - mean(parcelBlockAveragedChrom(iParc, iBlock, :), 'omitnan');
                else
                    params.blockRemoved = iBlock; %will either be first or last
                end
            end
        end
        
        parcelBlockAveraged{iChrom} = parcelBlockAveragedChrom;

%         figure; plot(squeeze(mean(parcelBlockAveraged(:, :, :), 1)))
%         for iParcel = 1:numParcels
%             figure;
%             plot(timeAxis, squeeze(mean(parcelBlockAveraged(iParcel, :, :), 2)));
%             xlabel(strcat('Parcel ', num2str(iParcel)));
%         end
    end
end