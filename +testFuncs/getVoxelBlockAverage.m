function [voxelBlockData, voxelParcelIDs, params] = getVoxelBlockAverage(parcelSensMask, chromData, info, params)

% Extract block-averaged, baseline-corrected data from every voxel

    % Get voxel locations in parcels
    [idx1, idx2, idx3] = ind2sub(size(parcelSensMask), find(parcelSensMask > 0));
    voxelParcelIDs = parcelSensMask(parcelSensMask > 0);
    numVoxels = length(voxelParcelIDs);
    numTimepoints = size(chromData, 4);

    % Time window (in samples)
    dtPre = floor(params.dtPre * (info.system.framerate / 10)); 
    dtAfter = floor(params.dtAfter * (info.system.framerate / 10));
    numSamples = dtPre + dtAfter + 1;

    % Get stimulus timepoints from Pulse_2–Pulse_6
    allPulses = [];
    for i = 2:6
        fieldName = sprintf('Pulse_%d', i);
        if isfield(info.paradigm, fieldName)
            allPulses = [allPulses; info.paradigm.(fieldName)(:)];
        end
    end
    allPulses = info.paradigm.synchpts(allPulses);
    [allPulses, params.blockRemoved] = analysisTools.chopBlocks(allPulses, numTimepoints, dtPre, dtAfter);
    numBlocks = length(allPulses);

    % Initialise output
    voxelBlockData = zeros(numVoxels, numBlocks, numSamples);

    % Loop through voxels
    for iVoxel = 1:numVoxels
        ts = squeeze(chromData(idx1(iVoxel), idx2(iVoxel), idx3(iVoxel), :));

        for iBlock = 1:numBlocks
            startIdx = allPulses(iBlock) - dtPre;
            endIdx = allPulses(iBlock) + dtAfter;

            if startIdx >= 1 && endIdx <= numTimepoints
                block = ts(startIdx:endIdx);
                voxelBlockData(iVoxel, iBlock, :) = block - mean(block, 'omitnan');
            else
                params.blockRemoved = [params.blockRemoved; iBlock];
            end
        end
    end
end

