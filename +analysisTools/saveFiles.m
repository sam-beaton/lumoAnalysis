function saveFiles(params, nirsFile, parcelData, cortexHb, cortexHbPeak, fooV, cortexMuA, channelData, tileData)

% Save outputs of nDotAnalysis.m

    [nirsFileDirec, nirsFileName, ~] = fileparts(nirsFile);
    dirParts = strsplit(nirsFileDirec, '/');
    nameParts = split(nirsFileName, '_');

    if exist('parcelData', 'var') && ~isempty(parcelData)
        % set filename
        parcelFileName = strcat(nameParts{1}, '_', nameParts{2}, '_', nameParts{3}, '_', nameParts{4}, '_', nameParts{5}, '_', 'parcelHb.mat');
        parcelDirectoryName = fullfile(params.outputDir, ...
            'parcelHb', ...
            dirParts{size(dirParts, 2)-2}, ...
            dirParts{size(dirParts, 2)-1}, ...
            dirParts{size(dirParts, 2)}, ...
            nameParts{4});
        parcelFileName = fullfile(parcelDirectoryName, parcelFileName);

        %save parcel data
        if ~isfolder(parcelDirectoryName)
            mkdir(parcelDirectoryName);
        end
        save(parcelFileName, 'parcelData');
    end

    if exist('cortexHb', 'var') && ~isempty(cortexHb)
        cortexHbFileName = strcat(nameParts{1}, '_', nameParts{2}, '_', nameParts{3}, '_', nameParts{4}, '_', 'cortexHb.mat');
        cortexHbDirectoryName = fullfile(params.outputDir, ...
            'cortexHb', ...
            dirParts{size(dirParts, 2)-2}, ...
            dirParts{size(dirParts, 2)-1}, ...
            dirParts{size(dirParts, 2)}, ...
            nameParts{4});
        cortexHbFileName = fullfile(cortexHbDirectoryName, cortexHbFileName);

        %save cortexHb
        if ~isfolder(cortexHbDirectoryName)
            mkdir(cortexHbDirectoryName);
        end
        save(cortexHbFileName, 'cortexHb');
    end

    if exist('cortexHbPeak', 'var') && ~isempty(cortexHbPeak)
        cortexHbPeakFileName = strcat(nameParts{1}, '_', nameParts{2}, '_', nameParts{3}, '_', nameParts{4}, '_', 'cortexHbPeak.mat');
        cortexHbPeakDirectoryName = fullfile(params.outputDir, ...
            'cortexHbPeak', ...
            dirParts{size(dirParts, 2)-2}, ...
            dirParts{size(dirParts, 2)-1}, ...
            dirParts{size(dirParts, 2)}, ...
            nameParts{4});
        cortexHbPeakFileName = fullfile(cortexHbPeakDirectoryName, cortexHbPeakFileName);

        %save cortexHb
        if ~isfolder(cortexHbPeakDirectoryName)
            mkdir(cortexHbPeakDirectoryName);
        end
        save(cortexHbPeakFileName, 'cortexHbPeak');
    end

    if exist('fooV', 'var') && ~isempty(fooV)
        fooVFileName = strcat(nameParts{1}, '_', nameParts{2}, '_', nameParts{3}, '_', nameParts{4}, '_', 'fooV.mat');
        fooVDirectoryName = fullfile(params.outputDir, ...
            'fooV', ...
            dirParts{size(dirParts, 2)-2}, ...
            dirParts{size(dirParts, 2)-1}, ...
            dirParts{size(dirParts, 2)}, ...
            nameParts{4});
        fooVFileName = fullfile(fooVDirectoryName, fooVFileName);

        %save cortexHb
        if ~isfolder(fooVDirectoryName)
            mkdir(fooVDirectoryName);
        end
        save(fooVFileName, 'fooV');
    end

    
    if exist('cortexMuA', 'var') && ~isempty(cortexMuA)
        cortexMuaFileName = strcat(nameParts{1}, '_', nameParts{2}, '_', nameParts{3}, '_', nameParts{4}, '_', 'cortexMuA.mat');
        cortexMuaDirectoryName = fullfile(params.outputDir, ...
            'cortexMuA', ...
            dirParts{size(dirParts, 2)-2}, ...
            dirParts{size(dirParts, 2)-1}, ...
            dirParts{size(dirParts, 2)}, ...
            nameParts{4});
        cortexMuaFileName = fullfile(cortexMuaDirectoryName, cortexMuaFileName);
        
        %save cortexMuA
        if ~isfolder(cortexMuaDirectoryName)
            mkdir(cortexMuaDirectoryName);
        end
        save(cortexMuaFileName, 'cortexMuA');
    end

    if exist('channelData', 'var') && ~isempty(channelData)
        % set filename
        channelFileName = strcat(nameParts{1}, '_', nameParts{2}, '_', nameParts{3}, '_', nameParts{4}, '_', nameParts{5}, '_', 'channelHb.mat');
        channelDirectoryName = fullfile(params.outputDir, ...
            'channelHb', ...
            dirParts{size(dirParts, 2)-2}, ...
            dirParts{size(dirParts, 2)-1}, ...
            dirParts{size(dirParts, 2)}, ...
            nameParts{4});
        channelFileName = fullfile(channelDirectoryName, channelFileName);

        %save parcel data
        if ~isfolder(channelDirectoryName)
            mkdir(channelDirectoryName);
        end
        save(channelFileName, 'channelData');
    end

    if exist('tileData', 'var') && ~isempty(tileData)
        % set filename
        tileFileName = strcat(nameParts{1}, '_', nameParts{2}, '_', nameParts{3}, '_', nameParts{4}, '_', nameParts{5}, '_', 'tileHb.mat');
        tileDirectoryName = fullfile(params.outputDir, ...
            'tileHb', ...
            dirParts{size(dirParts, 2)-2}, ...
            dirParts{size(dirParts, 2)-1}, ...
            dirParts{size(dirParts, 2)}, ...
            nameParts{4});
        tileFileName = fullfile(tileDirectoryName, tileFileName);

        %save parcel data
        if ~isfolder(tileDirectoryName)
            mkdir(tileDirectoryName);
        end
        save(tileFileName, 'tileData');
    end
end