function trialNumbers = getTrialNumbers(info)

% Returns trial numbers based on paradigmFull info for recording

    synchtypes = info.paradigmFull.synchtype;
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
    
%     % Check forward match is correct i.e. starts from beginning of HaND
%     if isequal(synchtypes, fullPattern(1:n))
%         trialNumbers = 1:n;
%         trialNumbers = trialNumbers';
%         return;
%     end
%     
%     % if not check reverse match (partial ending - second half of a HaND recording)
%     if n >= 10 && isequal(synchtypes(end-10:end), fullPattern(end-10:end))
%         trialNumbers = (25 - n + 1):25;
%         trialNumbers = trialNumbers';
%         return;
%     end
%     
%     % No valid match
%     warning('synchtypes does not match expected pattern. Skipping...\n');
%     trialNumbers = []; % Return empty 

end