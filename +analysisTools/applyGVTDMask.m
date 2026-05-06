function fcClean = applyGVTDMask(parcelAveraged, TkeepGVTD, info)
%APPLYGVTDMASK  Restrict parcel timeseries to GVTD-included samples within
%               the functional-connectivity task window.
%
%   fcClean = analysisTools.applyGVTDMask(parcelAveraged, TkeepGVTD, info)
%
% INPUTS
%   parcelAveraged : cell array (one entry per chromophore). Each cell is a
%                    [parcels x timepoints] matrix.
%   TkeepGVTD      : [timepoints x 1] vector of 1s (included) and 0s
%                    (excluded by GVTD scrubbing).
%   info           : struct with the fields:
%                       info.paradigmFull.synchpts    [1x2] start/end
%                                                      sample INDEXES of
%                                                      the FC task
%                       info.paradigmFull.synchtimes  [1x2] start/end
%                                                      TIMES (s) of the
%                                                      FC task
%                       info.system.framerate         sampling freq (Hz)
%
% OUTPUT
%   fcClean : struct with fields:
%       .data            cell array, same size as parcelAveraged. Each
%                        cell is [parcels x nKept] — only retained samples.
%       .keepMask        [timepoints x 1] logical. true = retained.
%       .blockIdx        [nKept x 1]. Block number (1,2,...) for each
%                        retained sample, identifying contiguous runs of
%                        included samples.
%       .blockDurations  [nBlocks x 1]. Duration of each block in seconds.
%       .adjStartIdx     First retained sample index (>= synchpts(1)).
%       .adjEndIdx       Last  retained sample index (<= synchpts(2)).
%       .adjStartTime    Time (s) corresponding to adjStartIdx.
%       .adjEndTime      Time (s) corresponding to adjEndIdx.
%       .framerate       Sampling frequency (Hz), passed through.
%
% Behaviour
%   1. The FC window is defined by info.paradigmFull.synchpts(1:2).
%      If the start sample is excluded by GVTD, the window start is
%      advanced to the first included sample after it. Similarly the
%      window end is retreated to the last included sample before it.
%   2. Samples outside the (adjusted) FC window are dropped. Inside the
%      window, only samples with TkeepGVTD == 1 are retained.
%   3. Each contiguous run of retained samples is given a block number,
%      and its duration in seconds is recorded — useful for dynamic
%      connectivity analyses where continuity matters.
%
% Throws an error if there are no GVTD-included samples inside the FC
% window.

    % ---- Pull out fields ----
    synchpts   = info.paradigmFull.synchpts;
    synchtimes = info.paradigmFull.synchtimes;
    fs         = info.system.framerate;

    tKeep = logical(TkeepGVTD(:));   % column vector
    nT    = numel(tKeep);

    % ---- Adjust window bounds to first/last included sample ----
    startIdx = synchpts(1);
    endIdx   = synchpts(2);

    adjStart = find(tKeep(startIdx:endIdx), 1, 'first');
    adjEnd   = find(tKeep(startIdx:endIdx), 1, 'last');

    if isempty(adjStart) || isempty(adjEnd)
        error('analysisTools:applyGVTDMask:NoIncludedSamples', ...
              'No GVTD-included samples in FC window.');
    end

    adjStartIdx = startIdx + adjStart - 1;
    adjEndIdx   = startIdx + adjEnd   - 1;

    % ---- Build final keep-mask ----
    windowMask                          = false(nT, 1);
    windowMask(adjStartIdx:adjEndIdx)   = true;
    keepMask                            = tKeep & windowMask;

    % ---- Identify contiguous included blocks ----
    blockIdxFull = zeros(nT, 1);
    d            = diff([false; keepMask; false]);
    blockStarts  = find(d ==  1);
    blockEnds    = find(d == -1) - 1;
    nBlocks      = numel(blockStarts);
    for iBlk = 1:nBlocks
        blockIdxFull(blockStarts(iBlk):blockEnds(iBlk)) = iBlk;
    end

    blockIdx       = blockIdxFull(keepMask);
    blockDurations = (blockEnds - blockStarts + 1) / fs;  % seconds

    % ---- Apply mask to each parcel timeseries ----
    % Each cell is parcels x timepoints; time is the 2nd dim.
    parcelAveragedKept = cell(size(parcelAveraged));
    for iChrom = 1:numel(parcelAveraged)
        parcelAveragedKept{iChrom} = parcelAveraged{iChrom}(:, keepMask);
    end

    % ---- Bundle outputs ----
    fcClean.data           = parcelAveragedKept;
    fcClean.adjStartIdx    = adjStartIdx;
    fcClean.adjEndIdx      = adjEndIdx;
    fcClean.adjStartTime   = (adjStartIdx - synchpts(1)) / fs + synchtimes(1);
    fcClean.adjEndTime     = (adjEndIdx   - synchpts(1)) / fs + synchtimes(1);
    fcClean.framerate      = fs;
    fcClean.keepMask       = keepMask;
    fcClean.blockIdx       = blockIdx;
    fcClean.blockDurations = blockDurations;
    
end