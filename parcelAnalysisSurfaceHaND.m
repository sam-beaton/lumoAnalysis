%% ndotWangFunctional.m
% Drop-in replacement for your existing parcellation script, using the
% Wang432 fine-grained functional parcellation instead of the coarser UNC
% volumetric parcellation.
%
% What changes vs your existing script:
%   - maskParc / infoParc loaded from Wang432_Reg_Head.nii.gz instead of
%     Parc_Reg_Head.nii.gz
%   - parcelData output is otherwise identical in structure — blockData,
%     blockAvgData, trialNumbers, parcelNumbers, capName all present
%   - parcelNumbers now range 1–430 (LH) and 431–860 (RH) instead of the
%     coarser UNC parcel IDs
%
% Everything downstream (getParcelSensitivity, getParcelAverageFull, block
% averaging, saving) is UNCHANGED.
%
% NOTE on hemisphere convention:
%   LH parcels: IDs  1 – 430
%   RH parcels: IDs 431 – 862
% This matches the offset applied in register_wang432.py.

%% Initialise and load toolboxes
clear all, close all;
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT-main'));
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer'));
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/lumoAnalysis'));

%% Pathing and parameters
% ---------- User-defined parameters ------------
timepoints = {'06', '12'}; % full: {'01', '06', '12'}
params.task      = 'hand';      % 'hand', 'fc1' or 'fc2'

driveName = '/Volumes/Extreme SSD/';

params.parentDir  = fullfile(driveName, 'dot');
params.preProcDir = 'imageRecon';
params.outputDir  = fullfile(params.parentDir, 'derivatives');

% Light and parcel sensitivity thresholds (unchanged)
params.lightSensitivityMin = 0.05;
params.parcPercentMin      = 0.5;

% ---------- Constant parameters  -------------
capNames   = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/capNames.csv';
capCSV     = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingData.csv';
jacobianDir = fullfile(driveName, 'imageRecon/neurodot/Jacobians/');

params.maxChannelDistance = 40;
params.dtPre              = 40;
params.dtAfter            = 190;
params.numSurfaceAvgTrials = 5;
params.peakTime           = 90;
params.peakDuration       = 40;

params.PD     = 1;
params.Scale  = 1;
params.Th.P   = 5e-2;
params.Th.N   = -params.Th.P;
params.TC     = 1;
params.cbmode = 0;
params.Cmap.P = colorcube;

% ---------- Derivative parameters -------------
params.dataLoc = fullfile(params.parentDir, 'derivatives', 'preproc');

% ---------- Load files needed for all iterations of the loop ----------
% % Age-specific head segmentation 
% if ~exist('maskSeg', 'var')
%     [maskSeg, infoSeg] = LoadVolumetricData( ...
%         strcat(params.timepoint, '_0Months3T_head_segVol'), ...
%         fullfile(driveName, strcat('imageRecon/neurodot/Segmentations/', ...
%                                    params.timepoint, 'mo')), ...
%         'nii');
% end

%% 
for iTime = 1:length(timepoints)

    params.timepoint = timepoints{iTime};
    clear prevCapName maskParc infoParc;

    % Start diary for this timepoint
    diaryFile = fullfile(params.outputDir, sprintf('parcelAnalysisSurface_%smo_%s.txt', ...
        params.timepoint, datestr(now, 'yyyymmdd_HHMMSS')));
    diary(diaryFile);
    diary on;
    fprintf('Diary started: %s\n', diaryFile);
    fprintf('\n========== Processing timepoint %s months ==========\n', params.timepoint);

    % Load Wang parcellation for this age
    % load Wang432 parcellation 
    % This file was produced by registerWang432.py and has the sam parent 
    % directory as existing Parc_Reg_Head.nii.gz files i.e. .../mri/registered
    [maskParc, infoParc] = LoadVolumetricData( ...
        strcat(params.timepoint, 'mo_WangFunctional_Reg_Head'), ...
        fullfile('/Volumes/Extreme SSD/', ...
                 strcat('mri/registered/UNC_Surf_to_NeuroDev/No Mask/', ...
                        params.timepoint, 'mo')), ...
        'nii.gz');
    lhParcels = unique(maskParc(maskParc >= 1   & maskParc <= 430));
    rhParcels = unique(maskParc(maskParc >= 431 & maskParc <= 862));
    fprintf('Wang mask loaded: %d LH parcels, %d RH parcels\n', numel(lhParcels), numel(rhParcels));

    % Search for task files
    matchingFiles = analysisTools.getAgeTaskNirsFiles(params);

    % Run analysis
    for nsub = 3:length(matchingFiles)
    
        [~, name, ~] = fileparts(matchingFiles{nsub});
        fprintf('\nAnalysing file %d: %s\n', nsub, name);
    
        paramsFile = params;
    
        try
    
            % ------ Load / get .nirs data -------
            data = analysisTools.getNdotFile(matchingFiles{nsub});
    
            % ---------- Get trial numbers ------------
            trialNumbers = analysisTools.getTrialNumbers(data.info);
            if isempty(trialNumbers)
                fprintf('synchtypes does not match expected pattern. Skipping...\n');
                continue;
            end
    
            % -------- Preprocessing ---------
            lmdata = logmean(data.d);
            paramsFile.keep = data.info.pairs.r3d < paramsFile.maxChannelDistance & ...
                              data.info.MEAS.GI;
            paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, ...
                                     paramsFile.timepoint, matchingFiles{nsub});
    
            data.info = FindGoodMeas(lmdata, data.info);
            lmdata    = detrend_tts(lmdata);
            [data.info.GVTD, data.info.DQ_metrics.med_GVTD] = ...
                CalcGVTD(lmdata(data.info.MEAS.GI & data.info.pairs.r2d < 20, :));
            gvtdThresh = data.info.DQ_metrics.med_GVTD + 5*mad(data.info.GVTD, 1);
            TkeepGVTD = data.info.GVTD <= gvtdThresh;
    
            % ------- Load Jacobian ------
            jacobFName = strcat('A_', paramsFile.capName, '_on_HD_Mesh_', ...
                                paramsFile.timepoint, 'mo.mat');
            if exist('prevCapName', 'var') && ...
                    numel(paramsFile.capName) == numel(prevCapName) && ...
                    all(paramsFile.capName == prevCapName)
                % reuse previous Jacobian
            else
                jacob = load(fullfile(jacobianDir, jacobFName), 'info', 'A');
            end
            jacob = analysisTools.reshapeJacob(jacob);
    
            % ------- Image reconstruction -------
            [cortexMuA, iJacob] = analysisTools.imageReconstruction( ...
                lmdata, jacob, data.info, paramsFile);
    
            % ------- Spectroscopy --------
            if ~exist('E', 'var')
                load('Extinction_Coefficients.mat');
                spectra{1} = 1; spectra{2} = 1;
                E = Generate_Spectroscopy_Matrix( ...
                    [data.info.pairs.lambda(1), data.info.pairs.lambda(end)], ...
                    spectra, ExtCoeffs.Prahl);
            end
            cortexHb = analysisTools.adaptedSpectroscopy_img(cortexMuA, E);
    
            % ------- Convert segmentation to DOT space -------
            % regMaskParc contains Wang432 IDs (1–860)
            % The call is identical to volumetric — only the content of maskParc has changed.
            regMaskParc = affine3d_img(maskParc, infoParc, jacob.info.tissue.dim, eye(4));
    
            % Dilate Wang parcellation labels to fill cortical volume.
            % The Wang parcellation is surface-derived so only labels a thin voxel
            % shell after resampling — dilation fills it out to match the volumetric
            % coverage of the coarse parcellation.
            % Strategy: dilate each label separately to avoid label bleeding across
            % parcel boundaries, then combine.
            fprintf('  Dilating Wang parcellation labels...\n');
            dilStart = tic;
            uniqueLabels = unique(regMaskParc(regMaskParc > 0));
            dilatedMask  = zeros(size(regMaskParc), 'int32');
            se           = strel('sphere', 3);
            
            for iLbl = 1:length(uniqueLabels)
                lbl     = uniqueLabels(iLbl);
                binMask = regMaskParc == lbl;
                dilBin  = imdilate(binMask, se);
                % Only write into voxels not yet claimed by an earlier label
                newVox  = dilBin & (dilatedMask == 0);
                dilatedMask(newVox) = lbl;
            end
            
            regMaskParc = double(dilatedMask);
            fprintf('  Dilation complete. Nonzero voxels: %d (%.1fs)\n', ...
                sum(regMaskParc(:) > 0), toc(dilStart));
    
            lhMask = unique(regMaskParc(regMaskParc >= 1   & regMaskParc <= 430));
            rhMask = unique(regMaskParc(regMaskParc >= 431 & regMaskParc <= 862));
            fprintf('Parcels in DOT space — LH: %d, RH: %d\n', numel(lhMask), numel(rhMask));
    
    %         % Add this temporarily after regMaskParc is computed
    %         figure;
    %         subplot(1,2,1)
    %         imagesc(max(regMaskParc >= 1 & regMaskParc <= 430, [], 3))
    %         title('LH parcels in DOT space (max projection)')
    %         colorbar
    %         
    %         subplot(1,2,2)
    %         imagesc(max(regMaskParc >= 431 & regMaskParc <= 862, [], 3))
    %         title('RH parcels in DOT space (max projection)')
    %         colorbar
    
            % ------- Light coverage -------
            fooV = analysisTools.getLightCoverage( ...
                jacob.A, iJacob, jacob.info, data.info, paramsFile);
    
            % After computing fooV
            sensMap = fooV.lambda1 > params.lightSensitivityMin;
            [sX, sY, sZ] = ind2sub(size(sensMap), find(sensMap));
            fprintf('Sensitivity map centroid: %.1f %.1f %.1f\n', mean(sX), mean(sY), mean(sZ));
            fprintf('Sensitivity X range: %.1f – %.1f\n', min(sX), max(sX));
            fprintf('Sensitivity Y range: %.1f – %.1f\n', min(sY), max(sY));
            fprintf('Sensitivity Z range: %.1f – %.1f\n', min(sZ), max(sZ));
            fprintf('Sensitivity map size: %d %d %d\n', size(sensMap));
            fprintf('Parcellation mask size: %d %d %d\n', size(regMaskParc));
    
    %         % Check spatial extent of sensitivity mask
    %         figure;
    %         subplot(1,2,1)
    %         imagesc(max(fooV.lambda1 > params.lightSensitivityMin & ...
    %             regMaskParc >= 1 & regMaskParc <= 430, [], 3))
    %         title('LH sensitive voxels')
    %         
    %         subplot(1,2,2)
    %         imagesc(max(fooV.lambda1 > params.lightSensitivityMin & ...
    %             regMaskParc >= 431 & regMaskParc <= 862, [], 3))
    %         title('RH sensitive voxels')
    
            % ------- Parcel sensitivity -------
            % getParcelSensitivity works on integer IDs — Wang432 IDs are handled
            % identically to the coarse UNC IDs.
            parcelsSens = analysisTools.getParcelSensitivity( ...
                regMaskParc, fooV, ...
                paramsFile.lightSensitivityMin, paramsFile.parcPercentMin);
    
            lhSens = unique(parcelsSens.lambda1(parcelsSens.lambda1 >= 1   & parcelsSens.lambda1 <= 430));
            rhSens = unique(parcelsSens.lambda1(parcelsSens.lambda1 >= 431 & parcelsSens.lambda1 <= 862));
            fprintf('Sensitive parcels — LH: %d, RH: %d\n', numel(lhSens), numel(rhSens));
    
            % ------- Convert Hb data to DOT volume --------
            regCortexHb      = cell(2, 1);
            regCortexHb{1}   = Good_Vox2vol(cortexHb(:, :, 1), jacob.info.tissue.dim);
            regCortexHb{2}   = Good_Vox2vol(cortexHb(:, :, 2), jacob.info.tissue.dim);
    
            % -------- Per-parcel averaging -----------
            % getParcelAverageFull is called identically — parcelNumbers in the
            % output will now be Wang432 IDs rather than coarse UNC IDs.
            parcelAveraged = cell(data.info.io.Nwl, 1);
            parcelNumbers  = cell(data.info.io.Nwl, 1);
            for iChrom = 1:numel(fieldnames(parcelsSens))
                lambdaField = ['lambda' num2str(iChrom)];
                [parcelAveraged{iChrom}, parcelNumbers{iChrom}] = ...
                    analysisTools.getParcelAverageFull( ...
                        parcelsSens.(lambdaField), regCortexHb{iChrom});
                parcelAveraged{iChrom} = detrend_tts(parcelAveraged{iChrom});
            end
    
            % ---------- Block averaging per parcel ----------
            parcelBlockData    = cell(data.info.io.Nwl, 1);
            parcelBlockAvgData = cell(data.info.io.Nwl, 1);
            allPulses = [];
            for i = 2:6
                fieldName = sprintf('Pulse_%d', i);
                if isfield(data.info.paradigm, fieldName)
                    thisPulse = data.info.paradigm.(fieldName);
                    allPulses = [allPulses; thisPulse(:)];
                end
            end

            allPulses = data.info.paradigm.synchpts(allPulses);
            blockLength = floor((data.info.system.framerate/10)*paramsFile.dtAfter);
            dtPre = floor((data.info.system.framerate/10)*paramsFile.dtPre);
            badTrials = false(length(allPulses), 1);
            for iTrial = 1:length(allPulses)
                windowStart = max(1, allPulses(iTrial) - dtPre);
                windowEnd = min(allPulses(iTrial) + blockLength - 1, length(TkeepGVTD));
                badTrials(iTrial) = any(~TkeepGVTD(windowStart:windowEnd));
            end
            if any(badTrials)
                fprintf('GVTD rejected %d / %d trials\n', sum(badTrials), length(allPulses));
                allPulses(badTrials) = [];
                trialNumbers(badTrials) = [];
            end
            fprintf('GVTD threshold: %.4f | Timepoints rejected: %d / %d (%.1f%%)\n', ...
                gvtdThresh, sum(~TkeepGVTD), length(TkeepGVTD), 100*mean(~TkeepGVTD));
            fprintf('Trials remaining after GVTD: %d / %d\n', length(allPulses), length(allPulses)+sum(badTrials));
    
            for iChrom = 1:numel(fieldnames(parcelsSens))
                [parcelBlockAvgData{iChrom}, ~, ~, parcelBlockData{iChrom}, blockRemoved] = ...
                    analysisTools.getDataAverageBlock( ...
                        parcelAveraged{iChrom}, ...
                        allPulses, ...
                        floor((data.info.system.framerate / 10) * paramsFile.dtPre), ...
                        floor((data.info.system.framerate / 10) * paramsFile.dtAfter));
            end
    
            paramsFile.blockRemoved = blockRemoved;
            if ~isempty(blockRemoved)
                trialNumbers(blockRemoved) = [];
            end
            fprintf('Blocks removed by chopBlocks: %d\n', length(blockRemoved));
            fprintf('Final trial count: %d\n', length(trialNumbers));
            fprintf('  DONE: %d LH sensitive parcels, %d RH sensitive parcels, %d trials\n', ...
                numel(lhSens), numel(rhSens), numel(trialNumbers));
    
            % ---------- Housekeeping ------------
            prevCapName = paramsFile.capName;

            % ---------- Get block data for each voxel ----------
            cortexBlockAvgData = cell(data.info.io.Nwl, 1);
            cortexHbPeak = cell(data.info.io.Nwl, 1);
            for iChrom = 1:numel(fieldnames(parcelsSens))
                [cortexBlockAvgData{iChrom}, ~, ~, cortexBlockData{iChrom}, ~] = ...
                    analysisTools.getDataAverageBlock( ...
                        cortexHb(:,:,iChrom), ...
                        allPulses, ...
                        floor((data.info.system.framerate/10)*paramsFile.dtPre), ...
                        floor((data.info.system.framerate/10)*paramsFile.dtAfter));
            end
            for iChrom = 1:numel(fieldnames(parcelsSens))
                cortexHbPeak{iChrom} = analysisTools.getCortexChromPeak(cortexBlockData{iChrom}, ...
                                                        data.info, ...
                                                        paramsFile, ...
                                                        trialNumbers);
            end
    
            % ---------- Join parcel data ----------
            parcelData                = struct;
            parcelData.blockData      = parcelBlockData;
            parcelData.blockAvgData   = parcelBlockAvgData;
            parcelData.trialNumbers   = trialNumbers;
            parcelData.parcelNumbers  = parcelNumbers;  % now Wang432 IDs
            parcelData.capName        = paramsFile.capName;
    
            % Hemisphere convenience fields — useful downstream
            % LH: parcel IDs 1–430, RH: 431–860
            parcelData.lhParcelIDs = parcelNumbers{1}(parcelNumbers{1} <= 430);
            parcelData.rhParcelIDs = parcelNumbers{1}(parcelNumbers{1} >  430) - 430;
    
            for iChrom = 1:2
                parcelData.blockData{iChrom} = permute( ...
                    parcelData.blockData{iChrom}, [1 3 2]);
            end
    
            % ---------- Save ----------
            % Note: if analysisTools.saveFiles uses a filename suffix based on the
            % parcellation, you may want to add a 'wang432' tag here to avoid
            % overwriting your existing coarse-parcel outputs.
            analysisTools.saveFiles(paramsFile, matchingFiles{nsub}, ...
                                [], parcelData, [], cortexHbPeak, fooV);
    
        catch ME
            fprintf('Could not run analysis for %s:\n  %s\n', name, ME.message);
        end
    
        clearvars -except params timepoints iTime jacob driveName capCSV capNames ...
                jacobianDir maskParc infoParc ...
                matchingFiles prevCapName E;
        close all;
    
    end % subject loop

    diary off;
    fprintf('Diary saved\n\n');

end % timepoint loop
