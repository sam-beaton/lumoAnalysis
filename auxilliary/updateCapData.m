redcapTable = readtable('/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingCommentsData.csv');
cappingTable = readtable('/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingData.csv');

redcapTable.redcap_repeat_instrument = [];
redcapTable.redcap_repeat_instance = [];

redcapTable.cap_position = strings(height(redcapTable), 1);

keyCols = {'con_participantid_q1', 'redcap_event_name', 'testing_id'};

% Find rows in redcapTable that match rows in cappingTable based on keyCols
[tf, loc] = ismember(redcapTable(:, keyCols), cappingTable(:, keyCols));

% Replace data for matched rows (replace ALL columns)
redcapTable(tf, :) = cappingTable(loc(tf), :);

%sort in ascending row order by testing_id
redcapTable = sortrows(redcapTable, 'testing_id');

writetable(redcapTable, '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingDataTEST.csv')