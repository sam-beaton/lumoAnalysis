function tKeep = scrubKeepBlocks(tKeep, params)

% Ensures all values from unwanted blocks/trials are removed from tKeep
% before image reconstruction

    tKeepRows = size(tKeep, 1);
    foundAllZeros = 0;
    
    while foundAllZeros == 0
    
        %initialise
        firstZeroIdx = zeros(tKeepRows, 1); 
        
        %find remainining excluded stim marker locations
        for iFind = 1:tKeepRows
            idx = find(tKeep(iFind,:) == 0, 1); % Find first zero
            if ~isempty(idx)
                firstZeroIdx(iFind) = idx;
            end
        end
    
        if sum(firstZeroIdx) == 0
            excludeIdx = tKeep == -101;
            tKeep(excludeIdx) = 0;
            foundAllZeros = 1;
        else
            % mark entire post-stimulus block
            for iFill = 1:tKeepRows
                if firstZeroIdx(iFill) ~=0
                    % fill block with temp variable for now so not identified on next loop
                    tKeep(iFill, firstZeroIdx(iFill):firstZeroIdx(iFill)+params.dtAfter) = -101; 
                end
            end
        end
    end

end