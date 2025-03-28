%% Pathing and parameters
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing - remove? only needed if generating meshes
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

%params
cohort = 'study';
timePoint = '01';
task = 'HAND';
preProcDir = 'cvPrune_NoSSR';
parentDir = '/Users/sambe/Documents/MATLAB/matlabProjects/PreP/data/indigo/';

filename = '3001_IN-C-016D_INDIGO_NIRS_1m_23-10-23_11-04-33_cvPrune_NoSSR.nirs';
capCSV = '/Users/sambe/Data_raw/indigo/study/studyLayoutFiles/cap_list'; %csv format
capInfo = readtable(capCSV);

testingIDs = capInfo.testing_id;
testingID = str2double(filename(1:4));
capName = cell2mat(capInfo{find(testingIDs == testingID), strcat('x', timePoint, 'mo_cap_name')});

jacobianDir=strcat('/Users/sambe/Data/imageRecon/neurodot/Jacobians/', capName, '_', timePoint, 'mo/');
meshDir = strcat('/Users/sambe/Data/imageRecon/neurodot/Meshes/', timePoint, 'mo/'); %folder containing meshes
outputDir=strcat('/Users/sambe/Data/imageRecon/neurodot/workbench/', capName, '_', timePoint, 'mo'); %Output Directory for files

nirsFile = fullfile(parentDir, cohort, strcat(timePoint, 'mo'), task, 'nirs', preProcDir, filename);
nDotFile = fullfile(parentDir, cohort, strcat(timePoint, 'mo'), task, 'nDot', preProcDir, filename);

%%
nirsFile = '/Users/sambe/dot/nirs/sub-053b/ses-12/nirs/sub-053b_ses-12_task-hand_run-02.nirs'
nDotFile = '/Users/sambe/dot/derivatives/ndot/sub-053b/ses-12/nirs/sub-053b_ses-12_task-hand_run-02'

%% Convert data to Neurodot format
%second variable = 1 saves file; change to 0 if saving not necessary
% data: (preprocessed) intensity/attenuation data
% info: SDfile/array file structure

% =============== NEED INFO IF LOADING MESH SPECIFIC LATER? =============
% ========== NEED FLAGS AS IN NEURODOT IMAGERECONSTRUCTION SCRIPT? ======

[data, info] = sbDataNirs2ndot(nirsFile, 1, nDotFile);

%% Load support files
% Load mesh
% Load parcellation 
% Load A-matrix
% Load PAD file
load([jacobianDir, 'Pad_HD_Mesh_', timePoint, 'mo_', capName, '.mat']); %loads as 'info'

%% View pre-processed data
keep = info.pairs.WL==2 & info.pairs.r2d < 40 & info.MEAS.GI; % measurements to include

figure('Position',[100 100 550 780])
subplot(3,1,1); plot(lmdata(keep,:)'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('log(\phi/\phi_0)') 
m=max(max(abs(lmdata(keep,:))));
subplot(3,1,2); imagesc(lmdata(keep,:),[-1,1].*m); 
colorbar('Location','northoutside');
xlabel('Time (samples)');ylabel('Measurement #')
[ftmag,ftdomain] = fft_tts(squeeze(mean(lmdata(keep,:),1)),info.system.framerate); % Generate average spectrum
subplot(3,1,3); semilogx(ftdomain,ftmag);
xlabel('Frequency (Hz)');ylabel('|X(f)|');xlim([1e-3 1])

nlrGrayPlots_180818(lmdata,info); % Gray Plot with synch points

%% Block Averaging the measurement data and view
badata = BlockAverage(lmdata, info.paradigm.synchpts(info.paradigm.Pulse_2), dt);

badata=bsxfun(@minus,badata,mean(badata,2));

figure('Position',[100 100 550 780])
subplot(2,1,1); plot(badata(keep,:)'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('log(\phi/\phi_0)') 
m=max(max(abs(badata(keep,:))));
subplot(2,1,2); imagesc(badata(keep,:),[-1,1].*m); 
colorbar('Location','northoutside');
xlabel('Time (samples)');ylabel('Measurement #')


%% RECONSTRUCTION PIPELINE
if ~exist('A', 'var')       % In case running by hand or re-running script
    A=load([strcat('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT/Support_Files/A_Matrices/',A_fn)],'info','A');
    if length(size(A.A))>2  % A data structure [wl X meas X vox]-->[(WL*meas) X vox]
        [Nwl,Nmeas,Nvox]=size(A.A);
        A.A=reshape(permute(A.A,[2,1,3]),Nwl*Nmeas,Nvox);
    end        
end
Nvox=size(A.A,2); %number of voxels in voxellated sensitivity map
Nt=size(lmdata,2); %number of time points in (resampled) recording
% initialise absorption matrix i.e. inversion of A, twice (once for each WL)
cortex_mu_a=zeros(Nvox,Nt,2); %(#voxels x #time points x #WLs)
for j = 1:2
    keep = (info.pairs.WL == j) & (info.pairs.r2d <= 40) & info.MEAS.GI; %find measurements for this wavelength less than 40mm in distance
    disp('> Inverting A')                
    iA = Tikhonov_invert_Amat(A.A(keep, :), 0.01, 0.1); % Invert A-Matrix
    disp('> Smoothing iA')
    iA = smooth_Amat(iA, A.info.tissue.dim, 3);         % Smooth Inverted A-Matrix      
    cortex_mu_a(:, :, j) = reconstruct_img(lmdata(keep, :), iA);% Reconstruct Image Volume
end

% Spectroscopy
%%%% ============ GENERATE VALUES FOR 'E' HERE USING sbPrePCalcDPF ========
if ~exist('E', 'var'),load('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT/Support_Files/A_Matrices/E.mat'),end
cortex_Hb = spectroscopy_img(cortex_mu_a, E);
cortex_HbO = cortex_Hb(:, :, 1);
cortex_HbR = cortex_Hb(:, :, 2);
cortex_HbT = cortex_HbO + cortex_HbR;

%% Select Volumetric visualizations of block averaged data
if ~exist('MNI', 'var')
    [MNI,infoB]=LoadVolumetricData('Segmented_MNI152nl_on_MNI111',[],'4dfp'); % load MRI (same data set as in A matrix dim)
end
MNI_dim = affine3d_img(MNI,infoB,A.info.tissue.dim,eye(4),'nearest'); % transform to DOT volume space 

%% Explore the block-averaged data interactively
badata_HbO = BlockAverage(cortex_HbO, info.paradigm.synchpts(info.paradigm.Pulse_2), dt); %block averaging of HbO data
badata_HbO=bsxfun(@minus,badata_HbO,badata_HbO(:,1)); %centres the data by subtrating the value at the first timepoint in each voxel
badata_HbOvol = Good_Vox2vol(badata_HbO,A.info.tissue.dim); %transforms B.A. data into a 4D volume

Params.Scale=0.8*max(abs(badata_HbOvol(:))); % used to scale colormapping
Params.Th.P=0.2*Params.Scale; % used to scale colormapping
Params.Th.N=-Params.Th.P; % used to scale colormapping
Params.Cmap='jet';
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,badata_HbOvol,info) %plot data in 4D (x,y,z,time)

%% Explore whole TS HbO data interactively
data_HbOvol = Good_Vox2vol(cortex_HbO,A.info.tissue.dim); %transforms entire data TS into a 4D volume
ParamsWhole = Params; %duplicate params from B.A.
ParamsWhole.Scale=0.8*max(abs(data_HbOvol(:))); %alter scale
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,ParamsWhole,data_HbOvol,info) %plot entire data TS in 4D (x,y,z,time)

%% Visualize block-averaged absorption data 
% Use the next line to select mu_a (wavelength of 1 or 2 in the third
% dimension)
badata_MuA = BlockAverage(cortex_mu_a(:,:,1), info.paradigm.synchpts(info.paradigm.Pulse_2), dt); %block averaging of absorption data
badata_MuA=bsxfun(@minus,badata_MuA,badata_MuA(:,1)); %centre the data
badata_MuAvol = Good_Vox2vol(badata_MuA,A.info.tissue.dim); %transforms B.A. data into a 4D volume

%Params for visualisation
ParamsMuA.Scale=0.8*max(abs(badata_MuAvol(:)));
ParamsMuA.Th.P=0.2*ParamsMuA.Scale;
ParamsMuA.Th.N=-ParamsMuA.Th.P;
ParamsMuA.Cmap='jet';

%Plot absorption data
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,badata_MuAvol,info)

%% Visualize block-averaged HbO, HbR, or HbT data 
% use next line to change to HbO, HbR, or HbT (first input argument
badata_Hb = BlockAverage(cortex_HbT, info.paradigm.synchpts(info.paradigm.Pulse_2), dt); %B.A. of haemoglobin data 
badata_Hb=bsxfun(@minus,badata_HbO,badata_Hb(:,1)); %centre data
badata_Hbvol = Good_Vox2vol(badata_Hb,A.info.tissue.dim); %transforms B.A. data into a 4D volume
%Params for visualisation
ParamsHb.Scale=0.8*max(abs(badata_HbOvol(:)));
ParamsHb.Th.P=0.2*ParamsHb.Scale;
ParamsHb.Th.N=-ParamsHb.Th.P;
ParamsHb.Cmap='jet';
% Plot Hb data
PlotSlicesTimeTrace(MNI_dim,A.info.tissue.dim,Params,badata_HbOvol,info)