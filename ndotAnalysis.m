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

% Run Analysis
for nsub = 1:length(matchingFiles) %01m: 59; 06mo: ? ; 12mo: 25

    [~, name, ~] = fileparts(matchingFiles{nsub});
    fprintf('\nAnalysing file %d: %s\n', nsub, name);

    % ------ Reset file-specific parameters -----
    paramsFile = params;

    try
    
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

        % measurements to include when plotting
        paramsFile.keep = data.info.pairs.r3d < paramsFile.maxChannelDistance & data.info.MEAS.GI; 
        % Get cap name 
        paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, paramsFile.timepoint, matchingFiles{nsub});

        data.info = FindGoodMeas(lmdata, data.info);                      % Detect Noisy Channels
        lmdata = detrend_tts(lmdata);                                             % Detrend Data
        [data.info.GVTD, data.info.DQ_metrics.med_GVTD] = CalcGVTD(lmdata(data.info.MEAS.GI & data.info.pairs.r2d<20,:));  
        
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
    
        % ------- Convert segmentation and parcellation to DOT volume space -------
        % convert segmentation 
        regMaskSeg = affine3d_img(maskSeg, infoSeg, jacob.info.tissue.dim, eye(4));
        % convert parcellation
        regMaskParc=affine3d_img(maskParc, infoParc, jacob.info.tissue.dim, eye(4)); 
         
        % ------- Find extent of array coverage for each wavelength -------
        % Calculate light coverage
        fooV = analysisTools.getLightCoverage(jacob.A, iJacob, jacob.info, data.info, paramsFile);
                
        % ------- Find intersection of array coverage and parcellation --------
        parcelsSens = analysisTools.getParcelSensitivity(regMaskParc, fooV, paramsFile.lightSensitivityMin, paramsFile.parcPercentMin);
        
        % ------- Convert HbO and HbR data to DOT volume space --------
        regCortexHb = cell(2,1);
        regCortexHb{1} = Good_Vox2vol(cortexHb(:, :, 1), jacob.info.tissue.dim); %HbO
        regCortexHb{2} = Good_Vox2vol(cortexHb(:, :, 2), jacob.info.tissue.dim); %HbR

        % -------- Take parcel average of chromophore data -----------
        parcelAveraged = cell(data.info.io.Nwl,1); %should be 2*1
        parcelNumbers = cell(data.info.io.Nwl,1);
        for iChrom = 1:numel(fieldnames(parcelsSens))
            lambdaField = ['lambda' num2str(iChrom)];
            [parcelAveraged{iChrom}, parcelNumbers{iChrom}] = analysisTools.getParcelAverageFull(parcelsSens.(lambdaField), regCortexHb{iChrom});
            parcelAveraged{iChrom} = detrend_tts(parcelAveraged{iChrom});
        end

        TkeepChrom = cell(data.info.io.Nwl,1);
        for iChrom = 1:numel(fieldnames(parcelsSens))
            [gvtd, gvtdMedian] = CalcGVTD(parcelAveraged{iChrom});
            tempTkeep = ones([size(gvtd)]);
            tempTkeep(gvtd > gvtdMedian + 5*mad(gvtd, 1)) = 0;
            TkeepChrom{iChrom} = tempTkeep;
        end
        
        % ---------- Get block data for each parcel ----------
        parcelBlockData = cell(data.info.io.Nwl,1);
        parcelBlockAvgData = cell(data.info.io.Nwl,1);
        allPulses = [];
        for i = 2:6
            fieldName = sprintf('Pulse_%d', i);
            if isfield(data.info.paradigm, fieldName)
                thisPulse = data.info.paradigm.(fieldName);
                allPulses = [allPulses; thisPulse(:)];
            end
        end
        allPulses = data.info.paradigm.synchpts(allPulses);
        for iChrom = 1:numel(fieldnames(parcelsSens))
            [parcelBlockAvgData{iChrom}, ~, ~, parcelBlockData{iChrom}] = testFuncs.getDataAverageBlock( ...
                                                    parcelAveraged{iChrom}, ...
                                                    allPulses, ...
                                                    floor((data.info.system.framerate/10)*paramsFile.dtPre), ...
                                                    floor((data.info.system.framerate/10)*paramsFile.dtAfter), ...
                                                    TkeepChrom{iChrom});
        end

        
        % ---------- Get block data for each voxel ----------
        TkeepChrom = cell(data.info.io.Nwl,1);
        for iChrom = 1:numel(fieldnames(parcelsSens))
            [gvtd, gvtdMedian] = CalcGVTD(cortexHb(:, :, 1));
            tempTkeep = ones([size(gvtd)]);
            tempTkeep(gvtd > gvtdMedian + 5*mad(gvtd, 1)) = 0;
            TkeepChrom{iChrom} = tempTkeep;
        end
        cortexBlockAvgData = cell(data.info.io.Nwl,1);
        cortexHbPeak = cell(data.info.io.Nwl,1);
        for iChrom = 1:numel(fieldnames(parcelsSens))
            [cortexBlockAvgData{iChrom}, ~, ~, cortexBlockData{iChrom}] = testFuncs.getDataAverageBlock( ...
                                                    cortexHb(:,:,iChrom), ...
                                                    allPulses, ...
                                                    floor((data.info.system.framerate/10)*paramsFile.dtPre), ...
                                                    floor((data.info.system.framerate/10)*paramsFile.dtAfter), ...
                                                    TkeepChrom{iChrom});
        end
        for iChrom = 1:numel(fieldnames(parcelsSens))
            cortexHbPeak{iChrom} = testFuncs.getCortexChromPeak(cortexBlockData{iChrom}, ...
                                                    data.info, ...
                                                    paramsFile, ...
                                                    trialNumbers);
        end
        
        % ---------- Checks and housekeeping ------------
        % for comparison with next cap to save loading Jacobian if possible
        prevCapName = paramsFile.capName; 
        
        % check all blocks were used for averaging i.e. weren't too close to 
        % ends of recording
        if isfield(paramsFile, 'blockRemoved')
            trialNumbers(paramsFile.blockRemoved) = []; %removes unused trial number (1st or last)
        end
        
        % ---------- Join parcel data ____________
        parcelData = struct;
        parcelData.blockData = parcelBlockData;
        parcelData.blockAvgData = parcelBlockAvgData;
        parcelData.trialNumbers = trialNumbers;
        parcelData.parcelNumbers = parcelNumbers;
        parcelData.capName = paramsFile.capName;
        % needed for Change point detection code in Python
        for iChrom = 1:2
            parcelData.blockData{iChrom} = permute(parcelData.blockData{iChrom}, [1 3 2]);
        end

        % ---------- Save variables ----------
        analysisTools.saveFiles(paramsFile, matchingFiles{nsub}, parcelData, [], cortexHbPeak, fooV);
        
    catch
        fprintf(strcat('Could not run analysis - look into manually.\n'))
    end

    % ---------- Tidying up ----------
    clearvars -except params jacob driveName capCSV capNames jacobianDir fooVDir maskSeg infoSeg maskParc infoParc matchingFiles prevCapName regCortexGroup cortexCount
    close all; %incase plotting used
end

%% TESTING


