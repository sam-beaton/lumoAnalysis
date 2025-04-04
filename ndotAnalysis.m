%% Initialise and load toolboxes
clear all, close all;
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT-main')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing - remove? only needed if generating meshes
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

%% Pathing and parameters
% ---------- User-defined parameters ------------
params.timepoint = '12'; %'01', '06' or '12'
params.task = 'hand'; %'hand', 'fc1' or 'fc2'

%overarching directory containing .nirs files
params.parentDir = '/Users/sambe/dot'; 
%processing method (%suffix after 'preproc-' in derivatives folder)
params.preProcDir = 'standard'; 
% directory for (statistical) outputs
outputDir=strcat('/Users/sambe/imageRecon/neurodot/workbench/'); %Output Directory for files
%cap info for each participant, in .csv format
capCSV = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingData.csv'; 

% ---------- Constant parameters -------------
% cap names corresponding to info in capCSV file
capNames = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/capNames.csv'; % DOESN'T CHANGE
% directory with Jacobians (A matrices) and PADs fitted to age-specific
% head models
jacobianDir=strcat('/Users/sambe/imageRecon/neurodot/Jacobians/');
% maximum channel distance to analyse
params.maxChannelDistance = 45; % Frijia et al (2021)
params.dtPre = 40; % block average start
params.dtAfter = 180; % block average finish

% ---------- Derivative parameters -------------
% mesh directory dependant on age only
meshDir = strcat('/Users/sambe/Data/imageRecon/neurodot/Meshes/', params.timepoint, 'mo/'); %folder containing meshes
%data location task and age dependant
params.dataLoc = fullfile(params.parentDir, 'derivatives', strcat('preproc-', params.preProcDir));

%% Search for task files 
matchingFiles = analysisTools.getAgeTaskNirsFiles(params);

% Run Analysis
for nsub = 25%1:length(matchingFiles) %01m: 59; 06mo: ? ; 12mo: 25

    % ------ Reset file-specific parameters -----
    paramsFile = params;
    
    % ------ load/get .nirs data in ndot file form -------
    [data, info] = analysisTools.getNdotFile(matchingFiles{nsub});
    %nirs = load(matchingFiles{nsub}, '-mat');
    %[data, info] = analysisTools.getNdotFile('/Users/sambe/dot/nirs/sub-053b/ses-12/nirs/sub-053b_ses-12_task-hand_run-02.nirs');

    % -------- Get data-dependent values and parameters ---------
    %for viewing preprocessed data
    lmdata = logmean(data.d); 
    % measurements to include when plotting
    paramsFile.keep = info.pairs.r3d < paramsFile.maxChannelDistance & info.MEAS.GI; 
    % Get cap name 
    paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, params.timepoint, matchingFiles{nsub});


    % TEST TEST TEST TEST TEST TEST
    %use lower lowpass cutoff
    %data.d = lowpass(data.d, 0.2, info.system.framerate);
    %detrend:
    %lmdata = detrend_tts(lmdata);
    %resample
    %[data.dcResampled, info] = analysisTools.adaptedResample_tts(data.dc, info, 1, 1e-5); % 1 Hz Resampling (1 Hz)
    
    % -------- View processed data -----------
    %analysisTools.viewProcessedData(lmdata, info, paramsFile);

    % ------- Calculate block averaged data ----------
    badata = analysisTools.adaptedBlockAverage(data.dc, params, info);

    % ------- View block averaged data ---------
    analysisTools.viewBlockAveraged(badata, paramsFile);
    
    % ------- Image reconstruction: Load Jacobian and invert for mu_a ------
    % get A matrix name
    jacobFName = strcat('A_', paramsFile.capName, '_on_HD_Mesh_', paramsFile.timepoint, 'mo.mat');

    % load A matrix
    if ~exist ('jacob', 'var')
        jacob=load(fullfile(jacobianDir, jacobFName),'info','A'); %load info and the matrix itself from the Jacobian
    end

    %reshape if necessary
    if length(size(jacob.A))>2  
        [nWL, nMeas, nVox]=size(jacob.A);
        jacob.A=reshape(permute(jacob.A,[2,1,3]),Nwl*nMeas,nVox);
    end

    nVox=size(jacob.A,2); %number of voxels in voxellated sensitivity map
    nT = size(data.dc,2); %number of time points in (resampled) recording
    nLambda = length(unique(info.pairs.lambda)); %number of wavelengths
    % initialise absorption matrix i.e. inversion of A, twice (once for each WL)
    cortexMuA = zeros(nVox, nT, nLambda); 

    % invert A
    for j = 1:nLambda
        keep = (info.pairs.WL == j) & (info.pairs.r3d <= params.maxChannelDistance) & info.MEAS.GI; 
        fprintf('Inverting Jacobian: wavelength %g\n', j)             
        iJacob = Tikhonov_invert_Amat(jacob.A(keep, :), 0.01, 0.1); % Invert A-Matrix
        fprintf('Smoothing Inverse: wavelength %g\n', j)   
        iJacob = analysisTools.adaptedSmoothAmat(iJacob, jacob.info.tissue.dim, 3); % Smooth Inverted A-Matrix
        fprintf('Reconstructing Volume: wavelength %g\n', j) 
        cortexMuA(:, :, j) = reconstruct_img(data.dc(keep, :), iJacob); % Reconstruct Image Volume
    end

    % ------- Spectroscopy --------
    % Generate extinction coeffts if not already in workspace
    if ~exist('E', 'var')
        %load('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT/Support_Files/A_Matrices/E.mat')
        %exs = GetExtinctions([nirs.SD.Lambda]);
        load('Extinction_Coefficients.mat'); % loads ext. coeffts, including Prahl's
        spectra{1}=1;spectra{2}=1; % placeholders
        %Generate matrix of ext. coeffts.
        % Rows: (1) lower lambda, (2) higher lambda. Cols: (1) HbO, (2) HbR
        E = Generate_Spectroscopy_Matrix([info.pairs.lambda(1), info.pairs.lambda(end)], spectra, ExtCoeffs.Prahl);
    end
    cortexHb = spectroscopy_img(cortexMuA, E);
    cortexHbO = cortexHb(:, :, 1);
    cortexHbR = cortexHb(:, :, 2);
    cortexHbT = cortexHbO + cortexHbR; %total
    cortexHbD = cortexHbO - cortexHbR; %difference

    % get parcellation name

    %clear jacob
end


%% Load support files
% Load mesh
% Load parcellation 
% Load A-matrix

%% stat testing
% close all;
% 
% keepSig = zeros(size(badataKeep, 1), 1);
% 
% for iChan = 1:size(badataKeep, 1)
%     [h, p] = ttest2(badataKeep(iChan, 1:(1+params.dtPre)), badataKeep(iChan, 100:140), 'Vartype', 'unequal'); % Welch's t-test
%     if h == 1
%         keepSig(iChan) = 1;
%     end
% end


% for iChan = 50:75 %size(badataKeep, 1)
%     figure
%     plot(badataKeep(iChan, :));
% end

%% Select Volumetric visualizations of block averaged data
if ~exist('MNI', 'var')
    [maskSeg,infoSeg]=LoadVolumetricData([strcat(params.timepoint,'_0Months3T_head_segVol')], ...
        strcat('/Users/sambe/imageRecon/neurodot/Segmentations/', params.timepoint,'mo'), ...
        'nii'); % load MRI (same data set as in A matrix dim)
end
t1 = affine3d_img(maskSeg, infoSeg, jacob.info.tissue.dim, eye(4)); % transform to DOT volume space 

%% Explore the block-averaged data interactively
badata_HbO = BlockAverage(cortexHbO, info.paradigmFull.synchpts(info.paradigmFull.Pulse_2), params.dtAfter); %block averaging of HbO data
badata_HbO=bsxfun(@minus,badata_HbO,badata_HbO(:,1)); %centres the data by subtrating the value at the first timepoint in each voxel
badata_HbOvol = Good_Vox2vol(badata_HbO,jacob.info.tissue.dim); %transforms B.A. data into a 4D volume

Params.Scale=0.8*max(abs(badata_HbOvol(:))); % used to scale colormapping
Params.Th.P=0.2*Params.Scale; % used to scale colormapping
Params.Th.N=-Params.Th.P; % used to scale colormapping
Params.Cmap='jet';
PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,Params,badata_HbOvol,info) %plot data in 4D (x,y,z,time)

%% Explore whole TS HbO data interactively
data_HbOvol = Good_Vox2vol(cortexHbO,jacob.info.tissue.dim); %transforms entire data TS into a 4D volume
ParamsWhole = Params; %duplicate params from B.A.
ParamsWhole.Scale=0.8*max(abs(data_HbOvol(:))); %alter scale
PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,ParamsWhole,data_HbOvol,info) %plot entire data TS in 4D (x,y,z,time)

%% Visualize block-averaged absorption data 
% Use the next line to select mu_a (wavelength of 1 or 2 in the third
% dimension)
badata_MuA = BlockAverage(cortexMuA(:,:,1), info.paradigm.synchpts(info.paradigm.Pulse_2), dt); %block averaging of absorption data
badata_MuA=bsxfun(@minus,badata_MuA,badata_MuA(:,1)); %centre the data
badata_MuAvol = Good_Vox2vol(badata_MuA,jacob.info.tissue.dim); %transforms B.A. data into a 4D volume

%Params for visualisation
ParamsMuA.Scale=0.8*max(abs(badata_MuAvol(:)));
ParamsMuA.Th.P=0.2*ParamsMuA.Scale;
ParamsMuA.Th.N=-ParamsMuA.Th.P;
ParamsMuA.Cmap='jet';

%Plot absorption data
PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,Params,badata_MuAvol,info)

%% Visualize block-averaged HbO, HbR, or HbT data 
% use next line to change to HbO, HbR, or HbT (first input argument
badata_Hb = BlockAverage(cortexHbT, info.paradigm.synchpts(info.paradigm.Pulse_2), dt); %B.A. of haemoglobin data 
badata_Hb=bsxfun(@minus,badata_HbO,badata_Hb(:,1)); %centre data
badata_Hbvol = Good_Vox2vol(badata_Hb,jacob.info.tissue.dim); %transforms B.A. data into a 4D volume
%Params for visualisation
ParamsHb.Scale=0.8*max(abs(badata_HbOvol(:)));
ParamsHb.Th.P=0.2*ParamsHb.Scale;
ParamsHb.Th.N=-ParamsHb.Th.P;
ParamsHb.Cmap='jet';
% Plot Hb data
PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,Params,badata_HbOvol,info)