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
params.preProcDir = 'standard'; 
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
params.maxChannelDistance = 45; % Frijia et al (2021): 45
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
        fullfile('/Users/sambe/TEMPSTORAGE/', strcat('mri/registered/UNC_to_NeuroDev/No Mask/', params.timepoint, 'mo')), ...
        'nii.gz');
end

%% Search for task files 
matchingFiles = analysisTools.getAgeTaskNirsFiles(params);

% Run Analysis
for nsub = 1:length(matchingFiles) %01m: 59; 06mo: ? ; 12mo: 25

    [~, name, ~] = fileparts(matchingFiles{nsub});
    fprintf('\nAnalysing file %d: %s\n', nsub, name);

    % ------ Reset file-specific parameters -----
    paramsFile = params;

    try
    
        % ------ load/get .nirs data in ndot file form -------
        data = analysisTools.getNdotFile(matchingFiles{nsub});
        
        % -------- Get data-dependent values and parameters ---------
        %for viewing preprocessed data & image recon/spectroscopy
        lmdata = logmean(data.d);
        % measurements to include when plotting
        paramsFile.keep = data.info.pairs.r3d < paramsFile.maxChannelDistance & data.info.MEAS.GI; 
        % Get cap name 
        paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, paramsFile.timepoint, matchingFiles{nsub});
        
        % -------- View processed data -----------
        %analysisTools.viewProcessedData(lmdata, data.info, paramsFile);
        
        % ------- Calculate block averaged data ----------
        %[badata, ~, ~, ~, tKeep] = analysisTools.adaptedBlockAverage(lmdata, paramsFile, data.info);
        
        % ------- View block averaged data ---------
        %analysisTools.viewBlockAveraged(badata, paramsFile);

        % ------- Find tKeep (useful if not running block average) --------
        tKeep = analysisTools.findtKeep(lmdata, paramsFile, data.info);
        
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
        %PlotSlices(t1, jacob.info.tissue.dim, paramsFile, fooV.lambda1) % HbO
        %PlotSlices(t1, jacob.info.tissue.dim, paramsFile, fooV.lambda2) % HbR
        
        % ------- Find intersection of array coverage and parcellation --------
        parcelsSens = analysisTools.getParcelSensitivity(regMaskParc, fooV, paramsFile.lightSensitivityMin, paramsFile.parcPercentMin);
        
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
        [parcelBlockData, paramsFile] = analysisTools.getParcelAverageBlock(parcelAveraged, data.info, paramsFile);
        
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
        parcelData.blockData = parcelBlockData;
        parcelData.trialNumbers = trialNumbers;
        parcelData.parcelNumbers = parcelNumbers;
        parcelData.capName = paramsFile.capName;
        
        % ---------- Save variables ----------
        analysisTools.saveFiles(paramsFile, matchingFiles{nsub}, parcelData);
        
    catch
        fprintf(strcat('Could not run analysis - look into manually.\n'))
    end

    % ---------- Tidying up ----------
    clearvars -except params jacob driveName capCSV capNames jacobianDir maskSeg infoSeg maskParc infoParc matchingFiles prevCapName
    close all; %incase plotting used
end

%% TESTING


