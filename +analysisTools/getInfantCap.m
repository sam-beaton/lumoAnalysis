function capName = getInfantCap(capCSV, capNames, timepoint, filename)

% Uses age and testing ID to obtain name of layout file of cap used during 
% testing
    
    % load cap info
    capInfoInfant = readtable(capCSV);
    capInfoNames = readtable(capNames);

    %use filename to find infant ID
    [~ , fileName, ~] = fileparts(filename);
    fileSplits = strsplit(fileName, '_');
    subSplit = strsplit(fileSplits{1}, '-'); %this is the 3 digit ID number
    
    % find matching IDs
    idMatchRows = contains(capInfoInfant.con_participantid_q1, subSplit{2}(1:3)); 
    % find correct age - underscores prevent 1/12 confusion
    idTable = capInfoInfant(idMatchRows, :);
    ageMatchRow = contains(idTable.redcap_event_name, strcat(timepoint(2), '_'));
    ageMatchRow = find(ageMatchRow);
    % get cap code usedin redcap
    capCode = num2str(idTable(ageMatchRow, :).cap_size_code);
    % use cap code to find layout filename
    capRow = find(capInfoNames.cap_code == str2num(capCode));
    capName = capInfoNames{capRow, "cap_name"}{1};

end