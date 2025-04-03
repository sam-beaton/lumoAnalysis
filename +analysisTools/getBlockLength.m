function [dtPre, dtAfter] = getBlockLength(info)

% Derive blocklength for block averaging from baseline stimuslus timing 
% info
   
    if isfield(info, 'paradigm') && isfield (info.paradigm, 'tHRF')
        dtAfter = abs(info.paradigm.tHRF(end)).*info.system.framerate;
        dtPre = abs(info.paradigm.tHRF(1)).*info.system.framerate;
    elseif isfield(info, 'paradigmFull') && isfield (info.paradigmFull, 'tHRF')
        dtAfter = abs(info.paradigmFull.tHRF(end)).*info.system.framerate;
        dtPre = abs(info.paradigmFull.tHRF(1)).*info.system.framerate;
    else
        % Get all baseline timings
        baseTimes = info.paradigmFull.synchtimes(info.paradigmFull.Pulse_1);
        baseTimes = baseTimes(2:end-1); %ignore end baseline stims in case of time delay
        
        % Initialise array to store differences
        baseDiffs = zeros(length(baseTimes)-1, 1);
    
        % Calculate time differences between baseline stim occurrences
        for iBase = 1:length(baseDiffs)-1
            baseDiffs(iBase) = baseTimes(iBase+1)-baseTimes(iBase);
        end
    
        % Take mean of differences and round down for post-stim block length
        dtAfter = floor(mean(baseDiffs));
    
        %calculate proportional pre-stim block for statistical significance
        %testing
        dtPre = ceil(dtAfter./5);
    
        % convert back to sample points scale
        dtAfter = dtAfter*info.system.framerate;
        dtPre = dtPre*info.system.framerate;
    end

end