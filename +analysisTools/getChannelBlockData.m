function [channelBlockData, sourceNumbers, detectorNumbers, params] = getChannelBlockData(dataIn, info, params)

% Generates and returns arrays of block channel data for both chromophores

    %% Parameters and Initialization.
    dims = size(dataIn); %'data': (channels)*(sample points)
    Nt = dims(end); % last sample point - Assumes time is always the last dimension.
    NDtf = (ndims(dataIn) > 2); % sets to one if data >2 dimensional (needs transforming)
    
    % account for different sampling frequencies
    dtPre = floor(params.dtPre*(info.system.framerate/10)); 
    dtAfter = floor(params.dtAfter*(info.system.framerate/10));
    
    %% N-D Input (for 3-D or N-D voxel spaces).
    if NDtf
        dataIn = reshape(dataIn, [], Nt);
    end
    
    % incorporate params.keep
    keepChans = 1:size(dataIn, 1);
    keepChans = keepChans(params.keep);
    
    srcNums = info.pairs.Src(params.keep);
    detNums = info.pairs.Det(params.keep);
    
    
    if isfield(info, 'paradigmFull') && isfield(info.paradigmFull, 'sCh') % should be, but check
    
        %find channel-specific negative stim. markers (motion)
        [samplepoints, channels] = find(info.paradigmFull.sCh < 0);
        negStims = [samplepoints, channels];
        negStims = sortrows(negStims);
        %restrict to non-baseline pulse/stimuli
        pulseRows = ismember(negStims(:,1), info.paradigmFull.synchpts(info.paradigmFull.synchtype ~= 1));
        negLocs = negStims(pulseRows, :);
        
        % find non-baseline pulse samples
        pulseSamples = info.paradigmFull.synchpts(info.paradigmFull.synchtype ~= 1);
    
        % check pulse not too close to start or end of recording 
        [pulseSamples, params.blockRemoved] = analysisTools.chopBlocks(pulseSamples, size(dataIn, 2), dtPre, dtAfter);
    
        % initialise outputs
        hbChanSplit = size(dataIn, 1)/2;
        hboChans = keepChans(keepChans <= hbChanSplit);
        hbrChans = keepChans(keepChans > hbChanSplit);
        
        hboData = zeros(length(hboChans), length(pulseSamples), (dtPre+dtAfter+1));
        hbrData = zeros(length(hbrChans), length(pulseSamples), (dtPre+dtAfter+1));
    
        for chan = 1:length(hboChans)   
            for pulseSample = 1:length(pulseSamples')
                if ~(any(all(negLocs == [pulseSamples(pulseSample), hboChans(chan)], 2)))
                    hboData(chan, pulseSample, :) = dataIn(hboChans(chan), pulseSamples(pulseSample)-dtPre:pulseSamples(pulseSample)+dtAfter);
                    hboData(chan, pulseSample, :) = hboData(chan, pulseSample, :) - hboData(chan, pulseSample, dtPre+1);
                else
                    hboData(chan, pulseSample, :) = NaN;
                end
            end
        end
        
        for chan = 1:length(hbrChans)   
            for pulseSample = 1:length(pulseSamples')
                if ~(any(all(negLocs == [pulseSamples(pulseSample), hbrChans(chan)], 2)))
                    hbrData(chan, pulseSample, :) = dataIn(hbrChans(chan), pulseSamples(pulseSample)-dtPre:pulseSamples(pulseSample)+dtAfter);
                    hbrData(chan, pulseSample, :) = hbrData(chan, pulseSample, :) - hbrData(chan, pulseSample, dtPre+1);
                else
                    hbrData(chan, pulseSample, :) = NaN;
                end
            end
        end
        
        % assign outputs
        %chrom data
        channelBlockData{1} = hboData;
        channelBlockData{2} = hbrData;
        %chrom-dependent optode numbers
        sourceNumbers{1} = srcNums(1:length(srcNums)/2);
        sourceNumbers{2} = srcNums(length(srcNums)/2+1:end);
        detectorNumbers{1} = detNums(1:length(detNums)/2);
        detectorNumbers{2} = detNums(length(detNums)/2+1:end);

        return
    
    else
        warning('sCh not contained in nDot file.')
        channelBlockData = [];
        sourceNumbers = [];
        detectorNumbers = [];
        return
    end

end


