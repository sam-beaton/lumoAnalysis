%% Initialise and load toolboxes
clear all, close all;
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT-main')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing - remove? only needed if generating meshes
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

%% Pathing and parameters
% ---------- User-defined parameters ------------
params.timepoint = '06'; %'01', '06' or '12'
params.task = 'hand'; %'hand', 'fc1' or 'fc2'

%storage drive - easier than changing all names all the time
driveName = '/Volumes/Extreme SSD/';

%overarching directory containing .nirs files
params.parentDir = fullfile(driveName, 'dot');
%processing method (%suffix after 'preproc-' in derivatives folder)
params.preProcDir = 'standard'; 
% directory for (statistical) outputs
params.outputDir = fullfile(params.parentDir, 'derivatives'); %Output Directory for files

%Light and parcel sensitivity thresholds for overlap
params.lightSensitivityMin = 0.05; %light
params.parcPercentMin = 0.5; % min. voxel coverage %age of parcels required, as a decimal

% maximum channel distance to analyse
params.maxChannelDistance = 45; % Frijia et al (2021): 45

% ---------- Constant parameters - no need to change -------------
% cap names corresponding to info in capCSV file
capNames = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/capNames.csv';
%cap info for each participant, in .csv format
capCSV = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingData.csv'; 
% directory with Jacobians (A matrices) and PADs fitted to age-specific
% head models
jacobianDir=fullfile(driveName, 'imageRecon/neurodot/Jacobians/');
% Block averaging
params.dtPre = 45; % start
params.dtAfter = 190; % finish
% Plotslices
params.PD=1;
params.Scale=1;
params.Th.P=5e-2;
params.Th.N=-params.Th.P;
params.TC = 1; %use true color mapping
params.cbmode = 0; %use custom colorbar limits
params.Cmap.P = colorcube; %use custom colormap

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
        fullfile(driveName, strcat('mri/registered/UNC_to_NeuroDev/No Mask/', params.timepoint, 'mo')), ...
        'nii.gz');
end

%% Search for task files 
matchingFiles = analysisTools.getAgeTaskNirsFiles(params);

% Run Analysis
for nsub = [1, 2, 3, 4, 6, 7]%1:length(matchingFiles) %01m: 59; 06mo: ? ; 12mo: 25

    [~, name, ~] = fileparts(matchingFiles{nsub});
    fprintf(strcat('\nAnalysing file: ', name, '\n'))

    % ------ Reset file-specific parameters -----
    paramsFile = params;

    try
    
        % ------ load/get .nirs data in ndot file form -------
        data = analysisTools.getNdotFile(matchingFiles{nsub});
        
        % -------- Get data-dependent values and parameters ---------
        %for viewing preprocessed data & image recon/spectroscopy
        lmdata = logmean(data.d);
        % define number of samples in block (defaults given above as a backup)
%         if isfield(data.info, 'paradigmFull') && isfield(data.info.paradigmFull, 'tHRF')
%             zeroLoc = find(data.info.paradigmFull.tHRF == 0);
%             paramsFile.dtPre = zeroLoc-1;
%             paramsFile.dtAfter = size(data.info.paradigmFull.tHRF,2)-paramsFile.dtPre;
%         end
        % measurements to include when plotting
        paramsFile.keep = data.info.pairs.r3d < paramsFile.maxChannelDistance & data.info.MEAS.GI; 
        % Get cap name 
        paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, paramsFile.timepoint, matchingFiles{nsub});
        
        % -------- View processed data -----------
        %analysisTools.viewProcessedData(lmdata, data.info, paramsFile);

        lmdata = lowpass(lmdata, params.lpCutoff, data.info.system.framerate);                     % Low Pass Filter (0.5 Hz)
        
        % ------- Calculate block averaged data ----------
        [badata, ~, ~, ~, tKeep] = analysisTools.adaptedBlockAverage(lmdata, paramsFile, data.info);
        
        % ------- View block averaged data ---------
        analysisTools.viewBlockAveraged(badata, paramsFile);
    catch
    end
end