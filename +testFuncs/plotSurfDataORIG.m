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
params.cortexDataDir = 'cortexHbPeak'; 
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
%params.numSurfaceAvgTrials = 5; % number of trials to average for surface images
params.peakTime = 90; %anticipated peak time for haemodynamic response
params.peakDuration = 40; %anticipated peak window duration
% PlotInterpSurfMesh
paramsPlot.CBar_on = 1; 
paramsPlot.PD = 0;


% ---------- Derivative parameters -------------
%data location task and age dependant
params.cortexDataLoc = fullfile(params.parentDir, 'derivatives', params.cortexDataDir);
%params.dataLoc = fullfile(params.parentDir, 'nirs');

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
matchingFiles = testFuncs.getAgeTaskCortexDataFiles(params);

%get count of files to be averaged
numParticipants = length(matchingFiles);

for nsub = 1:length(matchingFiles)

    [~, name, ~] = fileparts(matchingFiles{nsub});
    fprintf('\nAnalysing file %d: %s\n', nsub, name);

    paramsFile = params;

    load(matchingFiles{nsub});
    cortexHbo = cortexHbPeak{1}{2};
    cortexHbr = cortexHbPeak{2}{2};

%     cortexHbo(isnan(cortexHbo)) = 0;
%     cortexHbr(isnan(cortexHbr)) = 0;

    paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, paramsFile.timepoint, matchingFiles{nsub});

    % ------- Image reconstruction: Load Jacobian and invert for MuA ------
    % construct Jacobian filename, load, and reshape (if needed)
    jacobFName = strcat('A_', paramsFile.capName, '_on_HD_Mesh_', paramsFile.timepoint, 'mo.mat');
    % load info for Jacobian matrix
    infoJacob=load(fullfile(jacobianDir, jacobFName),'info'); %load info and the matrix itself from the Jacobian
    
%     fooVname = 
%     fooV = 
    % ---------- Transfrom FOV to space for mesh projection ---------
%     fooV.lambda1 = affine3d_img(fooV.lambda1, jacob.info.tissue.dim, infoSeg);
%     fooV.lambda2 = affine3d_img(fooV.lambda2, jacob.info.tissue.dim, infoSeg);

    hboVolData = Good_Vox2vol(cortexHbo, infoJacob.info.tissue.dim);
    hbrVolData = Good_Vox2vol(cortexHbr, infoJacob.info.tissue.dim);

    hboSurfData = affine3d_img(hboVolData, infoJacob.info.tissue.dim, infoSeg);
    hbrSurfData = affine3d_img(hbrVolData, infoJacob.info.tissue.dim, infoSeg);

    if ~exist('groupHboSurfData', 'var')
        % save original dimensions of surface data for reshaping
        surfD1 = size(hboSurfData, 1);
        surfD2 = size(hboSurfData, 2);
        surfD3 = size(hboSurfData, 3);

        %initialise arrays for t-tests and averaging
        groupHboSurfData = nan([surfD1*surfD2*surfD3 length(matchingFiles)]);
        groupHbrSurfData = nan([surfD1*surfD2*surfD3 length(matchingFiles)]);

        %fill voxel-wise data for this participant
        groupHboSurfData(:, nsub) = hboSurfData(:);
        groupHbrSurfData(:, nsub) = hbrSurfData(:);

    else
        %fill voxel-wise data for this participant
        groupHboSurfData(:, nsub) = hboSurfData(:);
        groupHbrSurfData(:, nsub) = hbrSurfData(:);
    end

%     if ~exist('groupHboSurfData', 'var')
%         groupHboSurfData = hboSurfData./numParticipants;
%         groupHbrSurfData = hbrSurfData./numParticipants;
%     else
%         groupHboSurfData = groupHboSurfData + hboSurfData./numParticipants;
%         groupHbrSurfData = groupHbrSurfData + hbrSurfData./numParticipants;
%     end
   
end

brainmeshLeft = load(strcat('/Users/sambe/TEMPSTORAGE/mri/meshes/leftHemisphereMesh', params.timepoint,'mo.mat')); % loads nodes and faces
brainmeshRight = load(strcat('/Users/sambe/TEMPSTORAGE/mri/meshes/rightHemisphereMesh', params.timepoint,'mo.mat'));

[hboH, hboP, hboCI, hboStats] = ttest(groupHboSurfData, 0, 'Dim', 2);
[hbrH, hbrP, hbrCI, hbrStats] = ttest(groupHbrSurfData, 0, 'Dim', 2);

hboStats.tstat = reshape(hboStats.tstat, surfD1, surfD2, surfD3);
hbrStats.tstat = reshape(hbrStats.tstat, surfD1, surfD2, surfD3);

groupHbtSurfData = groupHboSurfData + groupHbrSurfData;
groupHbdSurfData = groupHboSurfData - groupHbrSurfData;

[hbtH, hbtP, hbtCI, hbtStats] = ttest(groupHbtSurfData, 0, 'Dim', 2);
[hbdH, hbdP, hbdCI, hbdStats] = ttest(groupHbdSurfData, 0, 'Dim', 2);

hbtStats.tstat = reshape(hbtStats.tstat, surfD1, surfD2, surfD3);
hbdStats.tstat = reshape(hbdStats.tstat, surfD1, surfD2, surfD3);

max1 = max(abs(hboStats.tstat(:)));
max2 = max(abs(hbrStats.tstat(:)));
max3 = max(abs(hbtStats.tstat(:)));
max4 = max(abs(hbdStats.tstat(:)));
globalMax = max([max1, max2, max3, max4]);
paramsPlot.Scale = globalMax;
paramsPlot.Th.P = 0; 
paramsPlot.Th.N = -paramsPlot.Th.P; 

PlotInterpSurfMesh(hboStats.tstat, brainmeshLeft, brainmeshRight, infoSeg, paramsPlot);
PlotInterpSurfMesh(hbrStats.tstat, brainmeshLeft, brainmeshRight, infoSeg, paramsPlot);
PlotInterpSurfMesh(hbtStats.tstat, brainmeshLeft, brainmeshRight, infoSeg, paramsPlot);
PlotInterpSurfMesh(hbdStats.tstat, brainmeshLeft, brainmeshRight, infoSeg, paramsPlot);