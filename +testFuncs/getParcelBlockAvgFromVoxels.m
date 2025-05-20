function [parcelBlockAveraged, parcelNumbers] = getParcelBlockAvgFromVoxels(voxelBlockData, voxelParcelIDs)

% Average voxel block responses into parcel block responses

    parcelNumbers = unique(voxelParcelIDs);
    numParcels = length(parcelNumbers);
    [~, numBlocks, numSamples] = size(voxelBlockData);

    parcelBlockAveraged = zeros(numParcels, numBlocks, numSamples);

    for iParcel = 1:numParcels
        thisParcel = parcelNumbers(iParcel);
        idxVoxels = find(voxelParcelIDs == thisParcel);
        parcelBlockAveraged(iParcel, :, :) = mean(voxelBlockData(idxVoxels, :, :), 1, 'omitnan');
    end
end
