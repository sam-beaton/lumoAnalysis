function saveFiles(params, nirsFile, parcelData, cortexHb, cortexMuA)

% Save outputs of nDotAnalysis.m 

    [nirsFileDirec, nirsFileName, ~] = fileparts(nirsFile);
    dirParts = strsplit(nirsFileDirec, '/');
    nameParts = split(nirsFileName, '_');

    if exist('parcelData', 'var')
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


    if exist('cortexHb', 'var')
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
    
    if exist('cortexMuA', 'var')
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
end