%% Initialise and load toolboxes
clear all, close all;
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT-main')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing - remove? only needed if generating meshes
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/lumoAnalysis')); %contains edited functions where necessary for use in image recon

%% Pathing and parameters
% ---------- User-defined parameters ------------
params.timepoint = '01'; %'01', '06' or '12'
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
params.TC = 1; %use true color mapping
params.cbmode = 0; %use custom colorbar limits
params.Cmap.P = colorcube; %use custom co§lormap

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
        fullfile('/Users/sambe/TEMPSTORAGE/', strcat('mri/registered/UNC_to_NeuroDev/No Mask/', params.timepoint, 'mo')), ...
        'nii.gz');
end

%% Search for task files 
matchingFiles = analysisTools.getAgeTaskNirsFiles(params);


[~, name, ~] = fileparts(matchingFiles{1});

% ------ Reset file-specific parameters -----
paramsFile = params;

% ------ load/get .nirs data in ndot file form -------
data = testFuncs.getNdotFile(matchingFiles{nsub}, paramsFile);

% ---------- Get trial numbers ------------
trialNumbers = analysisTools.getTrialNumbers(data.info);
%check trial numbers make sense - continue to next file if not
if isempty(trialNumbers)
    fprintf('synchtypes does not match expected pattern. Skipping...\n');
    continue; % Skip to next iteration
end

% -------- Get data-dependent values and parameters ---------
%for viewing preprocessed data & image recon/spectroscopy
lmdata = logmean(data.d);

% load tables
activationTable = readtable('/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/Chapters_and_Writing/Thesis/8_FDA_INDiGO/tables/channel_t_test_results_corrected_SIGONLY.csv');
exthabTable = readtable('/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/Chapters_and_Writing/Thesis/8_FDA_INDiGO/tables/channel_t_test_exthab_corrected_SIGONLY.csv');

%split tables by chrom and age
%activation
act12Table = activationTable(activationTable.Age == 12, :);
act6Table = activationTable(activationTable.Age == 6, :);
act6HbOTable = act6Table(act6Table.Chrom == 1, :);
act6HbRTable = act6Table(act6Table.Chrom == 2, :);
act12HbOTable = act12Table(act12Table.Chrom == 1, :);
act12HbRTable = act12Table(act12Table.Chrom == 2, :);
%chromophore
extHabHbOTable = exthabTable(exthabTable.Chrom == 1, :);
extHabHbRTable = exthabTable(exthabTable.Chrom == 2, :);

%get source and detector columns
act6HbOSrc = act6Table.Source;
act6HbRSrc = act6Table.Detector;
act12HbOSrc = act12Table.Source;
act12HbRDet = act12Table.Detector;
extHabHbOSrc = extHabHbOTable.Source;
extHabHbODet = extHabHbOTable.Detector;
