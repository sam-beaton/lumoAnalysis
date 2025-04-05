function parcelsAveraged = getParcelAverageFull(parcelSensMask, chromData)

% Returns parcel averaged time series data across the whole recording

    % get unique parcels in sensitivity map
    parcelNumbers = unique(parcelSensMask);
    parcelNumbers(parcelNumbers == 0) = []; %disregard 0 values i.e. not a parcel
    
    %store number of samples
    numSamps = size(chromData, 4);

    % initialise parcelwise timeseries storage
    parcelsAveraged = zeros(length(parcelNumbers), numSamps);

    for iParcel = 1:length(parcelNumbers)
        %indices of current parcel
        [idx1, idx2, idx3] = ind2sub(size(parcelSensMask), find(parcelSensMask == parcelNumbers(iParcel)));

        numMatches = length(idx1);
        parcelData = zeros(numMatches, numSamps);

        for iTime = 1:numSamps
            for iMatch = 1:numMatches
                parcelData(iMatch, iTime) = chromData(idx1(iMatch), idx2(iMatch), idx3(iMatch), iTime);
            end
        end
        %store mean parcel data
        parcelsAveraged(iParcel, :) = mean(parcelData, 1);
    end
end