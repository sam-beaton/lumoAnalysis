function [pulseSamples, blockRemoved] = chopBlocks(pulseSamples, dataLength, dtPre, dtAfter)

% Removes pulses from list due to be saved/averaged and returns a list 

    % check pulse not too close to start or end of recording 
    % initialise tracking variables
    blockRemoved = [];
    originalIndices = 1:length(pulseSamples);
    
    % check all pulses have acceptable buffer to ends of file
    iRemoved = true;
    while iRemoved
        iRemoved = false;
        % check first pulse
        if pulseSamples(1) - dtPre < 1 || pulseSamples(1) + dtAfter > dataLength
            blockRemoved(end+1) = originalIndices(1);  % record original index
            pulseSamples(1) = [];  % remove pulse
            originalIndices(1) = [];  % update tracking
            iRemoved = true;
            continue;  % start again
        end
    
        % check last pulse
        if pulseSamples(end) - dtPre < 1 || pulseSamples(end) + dtAfter > dataLength
            blockRemoved(end+1) = originalIndices(end);  
            pulseSamples(end) = []; 
            originalIndices(end) = [];  
            iRemoved = true;
            % loop will recheck without continue
        end
    end
end