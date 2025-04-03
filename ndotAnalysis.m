%% Initialise and load toolboxes
clear all, close all;
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox
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

% ---------- Derivative parameters -------------
% mesh directory dependant on age only
meshDir = strcat('/Users/sambe/Data/imageRecon/neurodot/Meshes/', params.timepoint, 'mo/'); %folder containing meshes
%data location task and age dependant
params.dataLoc = fullfile(params.parentDir, 'derivatives', strcat('preproc-', params.preProcDir));

%% Search for task files 
matchingFiles = analysisTools.getAgeTaskNirsFiles(params);

% Run Analysis
for nsub = 5%25%1:length(matchingFiles)
    
    %load/get .nirs data in ndot file form
    [data, info] = analysisTools.getNdotFile(matchingFiles{nsub});
    %nirs = load(matchingFiles{nsub}, '-mat');
    %[data, info] = analysisTools.getNdotFile('/Users/sambe/dot/nirs/sub-053b/ses-12/nirs/sub-053b_ses-12_task-hand_run-02.nirs');

    lmdata = logmean(data); %for viewing preprocessed data
    keep = info.pairs.WL==2 & info.pairs.r2d < params.maxChannelDistance & info.MEAS.GI; % measurements to include


    %derive blocklength from stim info
    [params.dtPre, params.dtAfter] = analysisTools.getBlockLength(info);
    
    % Get cap name 
    capName = analysisTools.getInfantCap(capCSV, capNames, params.timepoint, matchingFiles{nsub});
    
    % get A matrix name
    A_fn = strcat('A_', capName, '_on_HD_Mesh_', params.timepoint, 'mo.mat');
    % load A matrix
    % invert A

    % get parcellation name


end


%% Load support files
% Load mesh
% Load parcellation 
% Load A-matrix


%% View pre-processed data
keep = info.pairs.WL==2 & info.pairs.r2d < params.maxChannelDistance & info.MEAS.GI; % measurements to include

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
close all;
badata = analysisTools.adaptedBlockAverage(lmdata, params, info);

badata=bsxfun(@minus,badata,mean(badata,2));

badataKeep = badata(keep,:);

figure('Position',[100 100 550 780])
subplot(2,1,1); plot(badataKeep'); 
set(gca,'XLimSpec','tight'), xlabel('Time (samples)'), 
ylabel('log(\phi/\phi_0)') 
m=max(max(abs(badata(keep,:))));
subplot(2,1,2); imagesc(badataKeep,[-1,1].*m); 
colorbar('Location','northoutside');
xlabel('Time (samples)');ylabel('Measurement #')

%% stat testing
close all;

keepSig = zeros(size(badataKeep, 1), 1);

for iChan = 1:size(badataKeep, 1)
    [h, p] = ttest2(badataKeep(iChan, 1:(1+params.dtPre)), badataKeep(iChan, 100:140), 'Vartype', 'unequal'); % Welch's t-test
    if h == 1
        keepSig(iChan) = 1;
    end
end


% for iChan = 50:75 %size(badataKeep, 1)
%     figure
%     plot(badataKeep(iChan, :));
% end



%% RECONSTRUCTION PIPELINE
if ~exist('A', 'var')       % In case running by hand or re-running script
    A=load(fullfile(jacobianDir, A_fn),'info','A');
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
    keep = (info.pairs.WL == j) & (info.pairs.r2d <= params.maxChannelDistance) & info.MEAS.GI; %find measurements for this wavelength less than 40mm in distance
    disp('> Inverting A')                
    iA = Tikhonov_invert_Amat(A.A(keep, :), 0.01, 0.1); % Invert A-Matrix
    disp('> Smoothing iA')
    iA = smooth_Amat(iA, A.info.tissue.dim, 3);         % Smooth Inverted A-Matrix      
    cortex_mu_a(:, :, j) = reconstruct_img(data(keep, :), iA);% Reconstruct Image Volume
end

%% Spectroscopy
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
cortex_Hb = spectroscopy_img(cortex_mu_a, E);
cortex_HbO = cortex_Hb(:, :, 1);
cortex_HbR = cortex_Hb(:, :, 2);
cortex_HbT = cortex_HbO + cortex_HbR;
cortex_HbD = cortex_HbO - cortex_HbR;

%% Select Volumetric visualizations of block averaged data
if ~exist('MNI', 'var')
    [maskSeg,infoSeg]=LoadVolumetricData([strcat(params.timepoint,'_0Months3T_head_segVol')], ...
        strcat('/Users/sambe/imageRecon/neurodot/Segmentations/', params.timepoint,'mo'), ...
        'nii'); % load MRI (same data set as in A matrix dim)
end
t1 = affine3d_img(maskSeg, infoSeg, A.info.tissue.dim, eye(4)); % transform to DOT volume space 

%% Explore the block-averaged data interactively
badata_HbO = BlockAverage(cortex_HbO, info.paradigmFull.synchpts(info.paradigmFull.Pulse_2), params.dtAfter); %block averaging of HbO data
badata_HbO=bsxfun(@minus,badata_HbO,badata_HbO(:,1)); %centres the data by subtrating the value at the first timepoint in each voxel
badata_HbOvol = Good_Vox2vol(badata_HbO,A.info.tissue.dim); %transforms B.A. data into a 4D volume

Params.Scale=0.8*max(abs(badata_HbOvol(:))); % used to scale colormapping
Params.Th.P=0.2*Params.Scale; % used to scale colormapping
Params.Th.N=-Params.Th.P; % used to scale colormapping
Params.Cmap='jet';
PlotSlicesTimeTrace(t1,A.info.tissue.dim,Params,badata_HbOvol,info) %plot data in 4D (x,y,z,time)

%% Explore whole TS HbO data interactively
data_HbOvol = Good_Vox2vol(cortex_HbO,A.info.tissue.dim); %transforms entire data TS into a 4D volume
ParamsWhole = Params; %duplicate params from B.A.
ParamsWhole.Scale=0.8*max(abs(data_HbOvol(:))); %alter scale
PlotSlicesTimeTrace(t1,A.info.tissue.dim,ParamsWhole,data_HbOvol,info) %plot entire data TS in 4D (x,y,z,time)

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
PlotSlicesTimeTrace(t1,A.info.tissue.dim,Params,badata_MuAvol,info)

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
PlotSlicesTimeTrace(t1,A.info.tissue.dim,Params,badata_HbOvol,info)