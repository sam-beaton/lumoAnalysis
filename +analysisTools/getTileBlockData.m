function [tileBlockData, tileNumbers] = getTileBlockData(channelBlockData, srcNums, detNums, params);
    
    % Extract capname (GA00XXX + optional _NF, remove orientation)
    tokens = regexp(params.capName, '^(GA00\d{3})(_NF)?(?:_(R|L|U|RU|LU))?$', 'tokens');
    tempCapName = [tokens{1}{1}, tokens{1}{2}];
    % obtain optodes in each tile
    [sourceTile, detTile] = analysisTools.getTileNumbers(tempCapName);

    % initiliase variables
    tileNumbers = [1:length(sourceTile)];
    hboData = zeros(length(sourceTile), size(channelBlockData{1}, 2), size(channelBlockData{1}, 3));
    hbrData = zeros(length(sourceTile), size(channelBlockData{1}, 2), size(channelBlockData{1}, 3));

    for tileIdx = 1:length(sourceTile)

        tileSources = sourceTile{tileIdx};
        tileDetectors = detTile{tileIdx};
        
        %HbO:
        %find rows (channels) where both source and detector are in the
        %tile (should be same, but calculating separately just in case)
        bothInTileHbO = ismember(srcNums{1,1}, tileSources) & ismember(detNums{1,1}, tileDetectors);
        bothInTileHbR = ismember(srcNums{1,2}, tileSources) & ismember(detNums{1,2}, tileDetectors);

        for blockNum = 1:size(channelBlockData{1}, 2)
            hboData(tileIdx, blockNum, :) = mean(channelBlockData{1}(bothInTileHbO, blockNum, :), 'omitnan');
            hbrData(tileIdx, blockNum, :) = mean(channelBlockData{2}(bothInTileHbR, blockNum, :), 'omitnan');
        end

    end

    %assign data to output var.
    tileBlockData{1} = hboData;
    tileBlockData{2} = hbrData;
end
