function tKeep = scrubKeepBlocks(tKeep, params, info)

% Ensures all values from unwanted blocks/trials are removed from tKeep
% before image reconstruction

    tKeepChannels = size(tKeep, 1);
    tKeepTimepoints = size(tKeep, 2);
    foundAllZeros = 0;

    % account for different sampling frequencies
    dtPre = floor(params.dtPre*(info.system.framerate/10)); 
    dtAfter = floor(params.dtAfter*(info.system.framerate/10));
    
    while foundAllZeros == 0
    
        %initialise
        firstZeroIdx = zeros(tKeepChannels, 1); 
        
        %find remainining excluded stim marker locations
        for iFind = 1:tKeepChannels
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
            for iFill = 1:tKeepChannels
                if firstZeroIdx(iFill) ~=0
                    % fill block with temp variable for now so not identified on next loop
                    tKeep(iFill, firstZeroIdx(iFill):min(firstZeroIdx(iFill)+dtAfter, tKeepTimepoints)) = -101; 
                end
            end
        end
    end

    tKeep(find(tKeep == -101)) = 0;

end