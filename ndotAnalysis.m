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

%storage drive - easier than changing all names all the time
driveName = '/Volumes/G-DRIVE ArmorATD/';
driveName = '/Users/sambe/';
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

%% Search for task files 
matchingFiles = analysisTools.getAgeTaskNirsFiles(params);

% Run Analysis
for nsub = 42%[7,27,28,40,49,56,66]%1:length(matchingFiles) %01m: 59; 06mo: ? ; 12mo: 25

    [~, name, ~] = fileparts(matchingFiles{nsub});
    fprintf(strcat('\nAnalysing file: ', name, '\n'))

    % ------ Reset file-specific parameters -----
    paramsFile = params;

    try
    
        % ------ load/get .nirs data in ndot file form -------
        data = analysisTools.getNdotFile(matchingFiles{nsub});
        %data = analysisTools.getNdotFile('/Volumes/G-DRIVE ArmorATD/dot/nirs/sub-087k/ses-01/nirs/sub-087k_ses-01_task-hand_run-01.nirs');
        
        % -------- Get data-dependent values and parameters ---------
        %for viewing preprocessed data & image recon/spectroscopy
        lmdata = logmean(data.d); 
        % measurements to include when plotting
        paramsFile.keep = data.info.pairs.r3d < paramsFile.maxChannelDistance & data.info.MEAS.GI; 
        % Get cap name 
        paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, params.timepoint, matchingFiles{nsub});
        
        % -------- View processed data -----------
        %analysisTools.viewProcessedData(lmdata, data.info, paramsFile);
        
        % ------- Calculate block averaged data ----------
        %[badata, ~, ~, ~, tKeep] = analysisTools.adaptedBlockAverage(lmdata, params, data.info);
        
        % ------- View block averaged data ---------
        %analysisTools.viewBlockAveraged(badata, paramsFile);

        % ------- Find tKeep (useful if not running block average) --------
        tKeep = analysisTools.findtKeep(lmdata, params, data.info);
        
        % ------- Image reconstruction: Load Jacobian and invert for MuA ------
        % construct Jacobian filename, load, and reshape (if needed)
        jacobFName = strcat('A_', paramsFile.capName, '_on_HD_Mesh_', paramsFile.timepoint, 'mo.mat');
        % load A matrix
        % if previous cap and position same as current one, same jacobian is used
        if exist('prevCapName', 'var') && ...
           numel(paramsFile.capName) == numel(prevCapName) && ...
           all(paramsFile.capName == prevCapName)
           %do nothing
        else 
            jacob=load(fullfile(jacobianDir, jacobFName),'info','A'); %load info and the matrix itself from the Jacobian
        end
        % Reshape if necessary
        jacob = analysisTools.reshapeJacob(jacob);
        
        % Remove all values from excluded blocks in the reconstruction data
        tKeep = analysisTools.scrubKeepBlocks(tKeep, paramsFile);
        lmdata = lmdata.*tKeep;
        
        % Perform image reconstruction
        [cortexMuA, iJacob] = analysisTools.imageReconstruction(lmdata, jacob, data.info, paramsFile);
        
        % ------- Spectroscopy --------
        % Generate extinction coeffts (E) if not already in workspace
        if ~exist('E', 'var')
            load('Extinction_Coefficients.mat'); % loads ext. coeffts, including Prahl's
            spectra{1}=1;spectra{2}=1; % placeholders
            % Rows: (1) lower lambda, (2) higher lambda. Cols: (1) HbO, (2) HbR
            E = Generate_Spectroscopy_Matrix([data.info.pairs.lambda(1), data.info.pairs.lambda(end)], spectra, ExtCoeffs.Prahl);
        end
        
        % Perform spectroscopy
        %fprintf('Performing Spectroscopy\n') 
        cortexHb = analysisTools.adaptedSpectroscopy_img(cortexMuA, E);
        
        %fprintf('Creating derivative chromophore matrices\n') 
        % Variables below needed? can just pass them as written.
        %cortexHbO = cortexHb(:, :, 1);
        %cortexHbR = cortexHb(:, :, 2);
        %cortexHbT = cortexHbO + cortexHbR; %total
        %cortexHbD = cortexHbO - cortexHbR; %difference
    
        % ------- Convert segmentation and parcellation to DOT volume space -------
        % convert segmentation 
        t1 = affine3d_img(maskSeg, infoSeg, jacob.info.tissue.dim, eye(4));
        % convert parcellation
        regMaskParc=affine3d_img(maskParc, infoParc, jacob.info.tissue.dim, eye(4)); 
        
        % ------- Find extent of array coverage for each wavelength -------
        % Calculate light coverage
        fooV = analysisTools.getLightCoverage(jacob.A, iJacob, jacob.info, data.info, paramsFile);
        % Plotting:
        %PlotSlices(t1, jacob.info.tissue.dim, params, fooV.lambda1) % HbO
        %PlotSlices(t1, jacob.info.tissue.dim, params, fooV.lambda2) % HbR
        
        % ------- Find intersection of array coverage and parcellation --------
        parcelsSens = analysisTools.getParcelSensitivity(regMaskParc, fooV, params.lightSensitivityMin, params.parcPercentMin);
        
        % Plotting:
        %PlotSlices(t1, jacob.info.tissue.dim, paramsFile, parcelsSens.lambda1)
        %PlotSlices(t1, jacob.info.tissue.dim, paramsFile, parcelsSens.lambda2)
        
        % ------- Convert HbO and HbR data to DOT volume space --------
        regCortexHb = cell(2,1);
        regCortexHb{1} = Good_Vox2vol(cortexHb(:, :, 1), jacob.info.tissue.dim); %HbO
        regCortexHb{2} = Good_Vox2vol(cortexHb(:, :, 2), jacob.info.tissue.dim); %HbR
        
        % -------- Take parcel average of chromophore data -----------
        parcelAveraged = cell(data.info.io.Nwl,1); %should be 2*1
        for iLambda = 1:numel(fieldnames(parcelsSens))
            lambdaField = ['lambda' num2str(iLambda)];
            [parcelAveraged{iLambda}, parcelNumbers] = analysisTools.getParcelAverageFull(parcelsSens.(lambdaField), regCortexHb{iLambda});
        end
        
        % ---------- Obtain block data for each parcel ----------
        parcelBlockAveraged = cell(data.info.io.Nwl,1); %should be 2*1
        for iLambda = 1:numel(fieldnames(parcelsSens))
            [parcelBlockAveraged{iLambda}, paramsFile] = analysisTools.getParcelAverageBlock(parcelAveraged, data.info, paramsFile);
        end
        
        % ---------- Get trial numbers ------------
        trialNumbers = analysisTools.getTrialNumbers(data.info);
        
        % ---------- Checks and housekeeping ------------
        % for comparison with next cap to save loading Jacobian if possible
        prevCapName = paramsFile.capName; 
    
        %check trial numbers make sense - continue to next file if not
        if isempty(trialNumbers)
            fprintf('synchtypes does not match expected pattern. Skipping...\n');
            continue; % Skip to next iteration
        end
        
        % check all blocks were used for averaging i.e. weren't too close to 
        % ends of recording
        if isfield(paramsFile, 'blockRemoved')
            trialNumbers(paramsFile.blockRemoved) = []; %removes unused trial number (1st or last)
        end
        
        % ---------- Join parcel data ____________
        parcelData = struct;
        parcelData.blockAverage = parcelBlockAveraged;
        parcelData.trialNumbers = trialNumbers;
        parcelData.parcelNumbers = parcelNumbers;
        parcelData.capName = paramsFile.capName;
        
        % ---------- Save variables ----------
        analysisTools.saveFiles(paramsFile, matchingFiles{nsub}, parcelData);
        
    catch
        fprintf(strcat('Could not run analysis - look into manually.\n'))
    end

    % ---------- Tidying up ----------
    clearvars -except params jacob driveName capCSV capNames jacobianDir meshDir maskSeg infoSeg maskParc infoParc matchingFiles prevCapName
    close all; %incase plotting used
end

%% TESTING


%% Load support files
% Load mesh



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


% %% Explore the block-averaged data interactively (HbO)
% badataHbO = BlockAverage(cortexHbO, data.info.paradigmFull.synchpts(data.info.paradigmFull.Pulse_2), params.dtAfter); %block averaging of HbO data
% badataHbO=bsxfun(@minus,badataHbO,badataHbO(:,1)); %centres the data by subtrating the value at the first timepoint in each voxel
% badataHbOvol = Good_Vox2vol(badataHbO,jacob.info.tissue.dim); %transforms B.A. data into a 4D volume
% 
% Params.Scale=0.8*max(abs(badataHbOvol(:))); % used to scale colormapping
% Params.Th.P=0.2*Params.Scale; % used to scale colormapping
% Params.Th.N=-Params.Th.P; % used to scale colormapping
% Params.Cmap='jet';
% PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,Params,badataHbOvol,data.info) %plot data in 4D (x,y,z,time)
% 
% %% Explore the block-averaged data interactively (HbD)
% badataHbD = BlockAverage(cortexHbD, data.info.paradigmFull.synchpts(data.info.paradigmFull.Pulse_2), params.dtAfter); %block averaging of HbD data
% badataHbD=bsxfun(@minus,badataHbD,badataHbD(:,1)); %centres the data by subtrating the value at the first timepoint in each voxel
% badataHbDvol = Good_Vox2vol(badataHbD,jacob.info.tissue.dim); %transforms B.A. data into a 4D volume
% 
% Params.Scale=0.8*max(abs(badataHbDvol(:))); % used to scale colormapping
% Params.Th.P=0.2*Params.Scale; % used to scale colormapping
% Params.Th.N=-Params.Th.P; % used to scale colormapping
% Params.Cmap='jet';
% PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,Params,badataHbDvol,data.info) %plot data in 4D (x,y,z,time)
% 
% %% Explore whole TS HbO data interactively
% data_HbOvol = Good_Vox2vol(cortexHbO,jacob.info.tissue.dim); %transforms entire data TS into a 4D volume
% ParamsWhole = Params; %duplicate params from B.A.
% ParamsWhole.Scale=0.8*max(abs(data_HbOvol(:))); %alter scale
% PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,ParamsWhole,data_HbOvol,data.info) %plot entire data TS in 4D (x,y,z,time)
% 
% %% Visualize block-averaged absorption data 
% % Use the next line to select mu_a (wavelength of 1 or 2 in the third
% % dimension)
% badata_MuA = BlockAverage(cortexMuA(:,:,1), data.info.paradigm.synchpts(data.info.paradigm.Pulse_2), dt); %block averaging of absorption data
% badata_MuA=bsxfun(@minus,badata_MuA,badata_MuA(:,1)); %centre the data
% badata_MuAvol = Good_Vox2vol(badata_MuA,jacob.info.tissue.dim); %transforms B.A. data into a 4D volume
% 
% %Params for visualisation
% ParamsMuA.Scale=0.8*max(abs(badata_MuAvol(:)));
% ParamsMuA.Th.P=0.2*ParamsMuA.Scale;
% ParamsMuA.Th.N=-ParamsMuA.Th.P;
% ParamsMuA.Cmap='jet';
% 
% %Plot absorption data
% PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,Params,badata_MuAvol,data.info)
% 
% %% Visualize block-averaged HbO, HbR, or HbT data 
% % use next line to change to HbO, HbR, or HbT (first input argument
% badata_Hb = BlockAverage(cortexHbT, data.info.paradigm.synchpts(data.info.paradigm.Pulse_2), dt); %B.A. of haemoglobin data 
% badata_Hb=bsxfun(@minus,badataHbO,badata_Hb(:,1)); %centre data
% badata_Hbvol = Good_Vox2vol(badata_Hb,jacob.info.tissue.dim); %transforms B.A. data into a 4D volume
% %Params for visualisation
% ParamsHb.Scale=0.8*max(abs(badataHbOvol(:)));
% ParamsHb.Th.P=0.2*ParamsHb.Scale;
% ParamsHb.Th.N=-ParamsHb.Th.P;
% ParamsHb.Cmap='jet';
% % Plot Hb data
% PlotSlicesTimeTrace(t1,jacob.info.tissue.dim,Params,badataHbOvol,data.info)