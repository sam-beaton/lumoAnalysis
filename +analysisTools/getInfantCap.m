function capName = getInfantCap(capCSV, capNames, timepoint, filename)

% Uses age and testing ID to obtain name of layout file of cap used during 
% testing

    try
    
        % load cap info
        capInfoInfant = readtable(capCSV); %cappingData.csv
        capInfoNames = readtable(capNames); %capNames.csv
    
        %use filename to find infant ID
        [~ , fileName, ~] = fileparts(filename);
        fileSplits = strsplit(fileName, '_');
        subSplit = strsplit(fileSplits{1}, '-'); %this is the 3 digit ID number
        
        % no frontal coverage for these IDs at 6mo - slightly different PAD files
        exceptionIDs = {'016', '018', '011', '002', '007', '020', '004', '001'}; 
    
        if ~ismember(subSplit{2}(1:3), exceptionIDs) || ~strcmp(timepoint, '06')
            % find matching IDs
            idMatchRows = contains(capInfoInfant.con_participantid_q1, subSplit{2}(1:3)); 
            % find correct age - underscores prevent 1/12 confusion
            idTable = capInfoInfant(idMatchRows, :);
            ageMatchRow = contains(idTable.redcap_event_name, strcat(timepoint(2), '_'));
            ageMatchRow = find(ageMatchRow);
            % get cap code usedin redcap
            capCode = num2str(idTable(ageMatchRow, :).cap_size_used);
            capPosition = idTable(ageMatchRow, :).cap_position{1};
            % use cap code to find layout filename
            capRow = find(capInfoNames.cap_code == str2num(capCode));
            capName = capInfoNames{capRow, "cap_name"}{1};
        
            if capPosition == 'S'
                capName = capInfoNames{capRow, "cap_name"}{1};
            elseif isempty(capPosition)
                capName = capInfoNames{capRow, "cap_name"}{1};
                warning('No cap position available: using default position.')
            else
                capName = capInfoNames{capRow, "cap_name"}{1};
                capName = strcat(capName, '_', capPosition);
            end
            return
        else
            switch subSplit{2}(1:3)
                case '016'
                    capName = 'GA00274';
                    return
                case {'018', '011', '002', '001'}
                    capName = 'GA00438_NF';
                    return
                case '020'
                    capName = 'GA00440_NF';
                    return
                case '007'
                    capName = 'GA00438_NF_L';
                    return
                case '004'
                    capName = 'GA00440_NF_U';
                    return
            end
        end
    catch
         warning('Issue fetching capName for this participant.');
         capName = []; % Return empty 
    end
end