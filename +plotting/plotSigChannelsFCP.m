%% Initialise and load toolboxes
clear all, close all;
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT-main')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing - remove? only needed if generating meshes
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/lumoAnalysis')); %contains edited functions where necessary for use in image recon

%% Pathing and parameters
% ---------- User-defined parameters ------------
params.timepoint = '12'; %'01', '06' or '12'
params.task = 'hand'; %'hand', 'fc1' or 'fc2'

%storage drive - easier than changing all names all the time
driveName = '/Volumes/Extreme SSD/';
    
%overarching directory containing .nirs files
params.parentDir = fullfile(driveName, 'dot');
%processing method (%suffix after 'preproc-' in derivatives folder)
params.preProcDir = 'imageRecon'; 
% directory for (statistical) outputs
params.outputDir = fullfile(params.parentDir, 'derivatives'); %Output Directory for files

%Light and parcel sensitivity thresholds for overlap
params.lightSensitivityMin = 0.05; %light
params.parcPercentMin = 0.5; % min. voxel coverage %age of parcels required, as a decimal

% ---------- Constant parameters - no need to change -------------
% cap names corresponding to info in capCSV file
capNames = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/capNames.csv';
%cap info for each participant, in .csv format
capCSV = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingData.csv'; 
% directory with Jacobians (A matrices) and PADs fitted to age-specific
% head models
jacobianDir=fullfile(driveName, 'imageRecon/neurodot/Jacobians/');
% maximum channel distance to analyse
params.maxChannelDistance = 40; % Frijia et al (2021): 45
% Block averaging
params.dtPre = 40; % start
params.dtAfter = 190; % finish
params.numSurfaceAvgTrials = 5; % number of trials to average for surface images
params.peakTime = 90; %anticipated peak time for haemodynamic response
params.peakDuration = 40; %anticipated peak window duration
% Plotslices
params.PD=1;
params.Scale=1;
params.Th.P=5e-2;
params.Th.N=-params.Th.P;
params.TC = 1; %use true colour mapping
params.cbmode = 0; %use custom colourbar limits
params.Cmap.P = colorcube; %use custom coulormap

% ---------- Derivative parameters -------------
%data location task and age dependant
params.dataLoc = fullfile(params.parentDir, 'derivatives', strcat('preproc-', params.preProcDir));

% ---------- Load files needed for all iterations of the loop ----------
% Age-specific head segmentation
if ~exist('maskSeg', 'var')
    [maskSeg,infoSeg]=LoadVolumetricData([strcat(params.timepoint,'_0Months3T_head_segVol')], ...
        fullfile(driveName, strcat('imageRecon/neurodot/Segmentations/', params.timepoint,'mo')), ...
        'nii');
end

% Age-specific cortical parcellation
if ~exist('maskParc', 'var')
    [maskParc,infoParc]=LoadVolumetricData([strcat(params.timepoint,'mo_Parc_Reg_Head')], ...
        fullfile('/Volumes/Extreme SSD/', strcat('mri/registered/UNC_to_NeuroDev/No Mask/', params.timepoint, 'mo')), ...
        'nii.gz');
end

%% Search for task files 
matchingFiles = analysisTools.getAgeTaskNirsFiles(params);


[~, name, ~] = fileparts(matchingFiles{1});

% ------ Reset file-specific parameters -----
paramsFile = params;

% ------ load/get .nirs data in ndot file form -------
data = testFuncs.getNdotFile(matchingFiles{1}, paramsFile);

% -------- Get data-dependent values and parameters ---------
%for viewing preprocessed data & image recon/spectroscopy
lmdata = logmean(data.d);

% load tables
activationTable = readtable('/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/Chapters_and_Writing/Thesis/8_FDA_INDiGO/tables/channel_valid_changepoints.csv');

% change chromophorecolumn name in table
activationTable.Properties.VariableNames{'Chromophore'} = 'Chrom';

%split tables by chrom and age
%activation
act12Table = activationTable(activationTable.Age == 12, :);
act6Table = activationTable(activationTable.Age == 6, :);
act6HbOTable = act6Table(act6Table.Chrom == 1, :);
act6HbRTable = act6Table(act6Table.Chrom == 2, :);
act12HbOTable = act12Table(act12Table.Chrom == 1, :);
act12HbRTable = act12Table(act12Table.Chrom == 2, :);


%get source and detector columns
act6HbOSrc = act6Table.Source;
act6HbRSrc = act6Table.Detector;
act12HbOSrc = act12Table.Source;
act12HbRDet = act12Table.Detector;


% set plotting params

colours.HbO  = [0.85, 0.33, 0.10];  % HbO = Red
colours.HbR  = [0, 0.45, 0.74];     % HbR = Blue
colours.Both = [0.47, 0.67, 0.19];  % Both = Green

styles.HbO   = '--';  % Red dashed
styles.HbR   = '-.';   % Blue dotted
styles.Both  = '-';   % Green solid

lineWidth = 2;  % Uniform line thickness for all activation lines

% ---------- plot 6mo actiuvation data

customPairs6.Src   = [];
customPairs6.Det   = [];
customPairs6.Color = {};
customPairs6.Style = {};
customPairs6.Width = [];

pairIDs_HbO = strcat(string(act6HbOTable.Source), '_', string(act6HbOTable.Detector));
pairIDs_HbR = strcat(string(act6HbRTable.Source), '_', string(act6HbRTable.Detector));
overlapPairs = intersect(pairIDs_HbO, pairIDs_HbR);

% HbO table
for i = 1:height(act6HbOTable)
    src = act6HbOTable.Source(i);
    det = act6HbOTable.Detector(i);
    pairID = strcat(string(src), '_', string(det));

    if ismember(pairID, overlapPairs)
        colour = colours.Both;
        style = styles.Both;
    else
        colour = colours.HbO;
        style = styles.HbO;
    end

    customPairs6.Src(end+1)   = src;
    customPairs6.Det(end+1)   = det;
    customPairs6.Color{end+1} = colour;
    customPairs6.Style{end+1} = style;
    customPairs6.Width(end+1) = lineWidth;
end

% HbR table
for i = 1:height(act6HbRTable)
    src = act6HbRTable.Source(i);
    det = act6HbRTable.Detector(i);
    pairID = strcat(string(src), '_', string(det));

    if ~ismember(pairID, pairIDs_HbO)
        colour = colours.HbR;
        style = styles.HbR;

        customPairs6.Src(end+1)   = src;
        customPairs6.Det(end+1)   = det;
        customPairs6.Color{end+1} = colour;
        customPairs6.Style{end+1} = style;
        customPairs6.Width(end+1) = lineWidth;
    end
end



% ----------- 12mo ativation plotting

customPairs12Act.Src   = [];
customPairs12Act.Det   = [];
customPairs12Act.Color = {};
customPairs12Act.Style = {};
customPairs12Act.Width = [];

% Create ID strings to detect overlap
pairIDs_HbO_12 = strcat(string(act12HbOTable.Source), '_', string(act12HbOTable.Detector));
pairIDs_HbR_12 = strcat(string(act12HbRTable.Source), '_', string(act12HbRTable.Detector));
overlapPairs_12 = intersect(pairIDs_HbO_12, pairIDs_HbR_12);

% Loop through HbO table
for i = 1:height(act12HbOTable)
    src = act12HbOTable.Source(i);
    det = act12HbOTable.Detector(i);
    pairID = strcat(string(src), '_', string(det));

    if ismember(pairID, overlapPairs_12)
        colour = colours.Both;
        style  = styles.Both;
    else
        colour = colours.HbO;
        style  = styles.HbO;
    end

    customPairs12Act.Src(end+1)   = src;
    customPairs12Act.Det(end+1)   = det;
    customPairs12Act.Color{end+1} = colour;
    customPairs12Act.Style{end+1} = style;
    customPairs12Act.Width(end+1) = lineWidth;
end

% Loop through HbR table
for i = 1:height(act12HbRTable)
    src = act12HbRTable.Source(i);
    det = act12HbRTable.Detector(i);
    pairID = strcat(string(src), '_', string(det));

    if ~ismember(pairID, pairIDs_HbO_12)
        colour = colours.HbR;
        style  = styles.HbR;

        customPairs12Act.Src(end+1)   = src;
        customPairs12Act.Det(end+1)   = det;
        customPairs12Act.Color{end+1} = colour;
        customPairs12Act.Style{end+1} = style;
        customPairs12Act.Width(end+1) = lineWidth;
    end
end


% ---------------- PLOTTING ------------

paramsPlot.dimension = '3D';

data.info.optodes.dpos2 = data.info.optodes.dpos3(:,1:2);
data.info.optodes.spos2 = data.info.optodes.spos3(:,1:2);

% Plot 6mo
testFuncs.testPlotCap(data.info, paramsPlot, customPairs6);

% Plot 12mo act
testFuncs.testPlotCap(data.info, paramsPlot, customPairs12Act);


