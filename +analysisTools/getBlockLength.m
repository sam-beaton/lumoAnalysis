function [dtPre, dtAfter] = getBlockLength(info)

% Derive blocklength for block averaging from baseline stimuslus timing 
% info
    
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
    dtPre = ceil(dt./10);

end