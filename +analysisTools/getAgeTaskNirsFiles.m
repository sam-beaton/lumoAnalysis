function matchingFiles = getAgeTaskNirsFiles(params)

% Fetches list of nirs files in given data location, returning them as a
% cell array

    % Recursively get all .nirs files with correct task and age
    timepointNum = str2double(params.timepoint(1:2)); % Converts '01' -> 1, '48' -> 48
    timepointNum = sprintf('%02d', timepointNum); % Ensures zero-padded format
    fileList = dir(fullfile(params.dataLoc, '**', sprintf('*ses-%s*task-%s*.nirs', timepointNum, params.task)));
    
    % Extract file paths directly
    matchingFiles = fullfile({fileList.folder}, {fileList.name});

end