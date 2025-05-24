filePath = '/Volumes/Extreme SSD/dot/derivatives/parcelHb';

fileList = dir(fullfile(filePath, '**', sprintf('*parcelHb.mat')));
fileList = fileList(arrayfun(@(f) ~startsWith(f.name, '.'), fileList));

% Extract file paths directly
matchingFiles = fullfile({fileList.folder}, {fileList.name});

for nsub = 1:length(matchingFiles)

    load(matchingFiles{nsub});

    if length(parcelData.trialNumbers) < 25

        parcelData.trialNumbers

    end

end


