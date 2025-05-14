rootDir = '/Volumes/Extreme SSD/dot/derivatives/parcelHb';
matFiles = dir(fullfile(rootDir, '**', '*.mat'));  % Recursively get all .mat files

% Convert to full paths
matFiles = matFiles(~startsWith({matFiles.name}, '.'));  % Omit hidden files like .DS_Store
matFilesFullPath = fullfile({matFiles.folder}, {matFiles.name});

for i = 1:length(matFilesFullPath)
    filePath = matFilesFullPath{i};
    load(filePath);  % Assumes variable 'parcelData' is loaded
    
    % Display which file is being loaded
    fprintf('Loaded file %d of %d: %s\n', i, length(matFilesFullPath), filePath);
    
    % Start a new, larger figure for each file
    fig = figure('Position', [100, 100, 1200, 600]);
    hold on;
    
    % Get number of parcels
    nParcels = size(parcelData.blockData{1}, 1);
    nTimepoints = size(parcelData.blockData{1}, 3);
    nBlocks = size(parcelData.blockData{1}, 2);
    
    for iParc = 1:nParcels
        % Get mean data for this parcel (averaged across blocks)
        parcelMean = squeeze(mean(parcelData.blockData{1}(iParc, :, :)));
        
        % Mean-center: subtract mean
        parcelMeanCentered = parcelMean - mean(parcelMean, 'omitnan');
        
        % Plot mean-centered data
        plot(parcelMeanCentered);
        
        % Now also update *all blocks* to be mean-centered
        for iBlock = 1:nBlocks
            originalData = parcelData.blockData{1}(iParc, iBlock, :);
            meanCenteredData = originalData - mean(originalData, 'omitnan');
            parcelData.blockData{1}(iParc, iBlock, :) = meanCenteredData;
        end
    end
    
    hold off;
    title(sprintf('File %d: %s', i, matFiles(i).name), 'Interpreter', 'none');
    xlabel('Time points');
    ylabel('Mean-centered Signal');
    legend(arrayfun(@(x) sprintf('Parcel %d', x), 1:nParcels, 'UniformOutput', false));
    
    % Close the figure
    close(fig);
    
    % Save back to the same .mat file (overwrites the file)
    save(filePath, 'parcelData');  % Use -v7.3 if data is large
    fprintf('Saved mean-centered data back to: %s\n', filePath);
end
