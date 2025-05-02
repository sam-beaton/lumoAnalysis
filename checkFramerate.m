%% Initialise and load toolboxes
clear all, close all;
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT-main')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing - remove? only needed if generating meshes
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

%% Pathing and parameters
% ---------- User-defined parameters ------------
params.timepoint = '01'; %'01', '06' or '12'
params.task = 'hand'; %'hand', 'fc1' or 'fc2'

%storage drive - easier than changing all names all the time
driveName = '/Volumes/Extreme SSD/';

%overarching directory containing .nirs files
params.parentDir = fullfile(driveName, 'dot');
%processing method (%suffix after 'preproc-' in derivatives folder)
params.preProcDir = 'standard'; 
% directory for (statistical) outputs
params.outputDir = fullfile(params.parentDir, 'derivatives'); %Output Directory for files
%cap info for each participant, in .csv format
capCSV = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingData.csv'; 

%Light and parcel sensitivity thresholds for overlap
params.lightSensitivityMin = 0.05; %light
params.parcPercentMin = 0.5; % min. voxel coverage %age of parcels required, as a decimal

% ---------- Constant parameters -------------
% cap names corresponding to info in capCSV file
capNames = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/capNames.csv'; % DOESN'T CHANGE
% directory with Jacobians (A matrices) and PADs fitted to age-specific
% head models
jacobianDir=fullfile(driveName, 'imageRecon/neurodot/Jacobians/');
% maximum channel distance to analyse
params.maxChannelDistance = 45; % Frijia et al (2021)
% Block averaging
params.dtPre = 40; % start
params.dtAfter = 180; % finish
% Plotslices
params.PD=1;
params.Scale=1;
params.Th.P=5e-2;
params.Th.N=-params.Th.P;
params.TC = 1; %use true color mapping
params.cbmode = 0; %use custom colorbar limits
params.Cmap.P = colorcube; %use custom colormap

% ---------- Derivative parameters -------------
% mesh directory dependant on age only
meshDir = fullfile(driveName, 'imageRecon/neurodot/Meshes/'); % ############# needed? #############
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

%%
ages = [6, 12];
for ageIdx = 1:length(ages)

    age = sprintf('%02d', ages(ageIdx));

    fprintf('CHECKING FRAMERATE FOR AGE %s\n', age);
    params.timepoint = age; %'01', '06' or '12'

    matchingFiles = analysisTools.getAgeTaskNirsFiles(params);
    
    % Run Analysis
    for nsub = 11:20%1:length(matchingFiles) %01m: 59; 06mo: ? ; 12mo: 25
    
        [~, name, ~] = fileparts(matchingFiles{nsub});
        fprintf(strcat('\nChecking file: ', name, '\n'))
    
        % ------ Reset file-specific parameters -----
        paramsFile = params;
    
        try
            
            nirs = load(matchingFiles{nsub}, '-mat');
            [row, col] = find(nirs.sCh == -2);
            A = [row, col];
            A = sortrows(A);
            framerate  = nirs.t(2)-nirs.t(1);
            fprintf(num2str(framerate))

            figure; plot(squeeze(nirs.dcAvg(:,1,:,2)))
            
        catch
            fprintf(strcat('Could not run analysis - look into manually.\n'))
        end
    
    end
end

%% TESTING


