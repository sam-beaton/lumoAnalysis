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
        % initialise tracking variables
        params.blockRemoved = [];
        originalIndices = 1:length(pulseSamples);
        
        % check all pulses have acceptable buffer to ends of file
        iRemoved = true;
        while iRemoved
            iRemoved = false;
            % check first pulse
            if pulseSamples(1) - dtPre < 1 || pulseSamples(1) + dtAfter > size(dataIn, 2)
                params.blockRemoved(end+1) = originalIndices(1);  % record original index
                pulseSamples(1) = [];  % remove pulse
                originalIndices(1) = [];  % update tracking
                iRemoved = true;
                continue;  % start again
            end
        
            % check last pulse
            if pulseSamples(end) - dtPre < 1 || pulseSamples(end) + dtAfter > size(dataIn, 2)
                params.blockRemoved(end+1) = originalIndices(end);  
                pulseSamples(end) = []; 
                originalIndices(end) = [];  
                iRemoved = true;
                % loop will recheck without continue
            end
        end
    
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
                else
                    hboData(chan, pulseSample, :) = -101;
                end
            end
        end
        
        for chan = 1:length(hbrChans)   
            for pulseSample = 1:length(pulseSamples')
                if ~(any(all(negLocs == [pulseSamples(pulseSample), hbrChans(chan)], 2)))
                    hbrData(chan, pulseSample, :) = dataIn(hbrChans(chan), pulseSamples(pulseSample)-dtPre:pulseSamples(pulseSample)+dtAfter);
                else
                    hbrData(chan, pulseSample, :) = -101;
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


