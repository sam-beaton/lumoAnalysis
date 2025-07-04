function data = getNdotFile(nirsFile, params)

% Uses a nirs filename (in BIDS format) to either load an existing
% corresponding ndot file if it exists, or create one if not

    %get component parts of filename to use for ndot file name
    [filePath, fileName, ~] = fileparts(nirsFile);
    derivsPath = regexp(filePath, '.*derivatives', 'match', 'once');
    derivsPath = fullfile(params.parentDir, 'derivatives');
    fNameSplits = strsplit(fileName, '_');
    taskSplit = strsplit(fNameSplits{3}, '-');
    %name nDot file
    nDotFile = fullfile(derivsPath, 'TESTndot', fNameSplits{1}, fNameSplits{2}, taskSplit{2}, fileName);
    
    %check existence of ndot file: load if exists, create if not
    if exist(strcat(nDotFile, '.mat'), 'file') == 2
        load(nDotFile)
    else
        % convert nirs to ndot file
        data = testFuncs.adaptedNirs2ndot(nirsFile, 1, nDotFile);
    end
end