function trialNumbers = getTrialNumbers(info)

% Returns trial numbers based on paradigmFull/paradigm info for recording
    
    if isfield(info, 'paradigmFull')
        synchtypes = info.paradigmFull.synchtype;
    else
        synchtypes = info.paradigm.synchtype;
    end
    synchtypes = synchtypes(synchtypes ~= 1)';
    % Define expected pattern
    fullPattern = repelem(2:6, 5); % [2,2,2,2,2,3,3,3,3,3,...,6,6,6,6,6]
    %find index where pattern match begins
    startID = strfind(fullPattern, synchtypes);
    
    if ~isempty(startID)
        trialNumbers = [startID:length(synchtypes) + startID - 1]';
        return
    else
        warning('synchtypes does not match expected pattern. Skipping...\n');
        trialNumbers = []; % Return empty
    end

end