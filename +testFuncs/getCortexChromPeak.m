function [cortexChromPeaks] = getCortexChromPeak(cortexChromBlocksData, info, params)

    % generate values for each voxel in the cortex, representative of the mean
    % of the peak of the first 'n' trials, where 'n' is a variable input
    
    peakRange = [params.peakTime - floor(params.peakDuration/2), params.peakTime + floor(params.peakDuration/2)];
    peakRange = floor(peakRange*(info.system.framerate/10));

    % take peak range from all voxels, trials
    dataPeakWindow = cortexChromBlocksData(:, peakRange(1):peakRange(2), 1:params.numSurfaceAvgTrials);

    % take mean of peak range values, per row
    trialPeaks = mean(dataPeakWindow, 2, 'omitnan');  % [voxels x 1 x trials]

    % take mean over trials (dim 3)
    cortexChromPeaks = mean(trialPeaks, 3, 'omitnan');  % [voxels x 1]

end