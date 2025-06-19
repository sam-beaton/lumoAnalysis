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
params.nBlocks = 6; %number of trial/condition blocks expected (including baseline)
%storage drive - easier than changing all names all the time
driveName = '/Volumes/Extreme SSD/';
    
%overarching directory containing .nirs files
params.parentDir = fullfile(driveName, 'dot');
%processing method (%suffix after 'preproc-' in derivatives folder)
params.cortexDataDir = 'cortexHbPeak';
params.parcelDataDir = 'parcelHb';
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
params.parcelDataLoc = fullfile(params.parentDir, 'derivatives', params.parcelDataDir);


% ---------- Load files needed for all iterations of the loop ----------
% Age-specific head segmentation
% if ~exist('maskSeg', 'var')
%     [maskSeg,infoSeg]=LoadVolumetricData([strcat(params.timepoint,'_0Months3T_head_segVol')], ...
%         fullfile(driveName, strcat('imageRecon/neurodot/Segmentations/', params.timepoint,'mo')), ...
%         'nii');
% end

% Age-specific cortical parcellation
if ~exist('maskParc', 'var')
    [maskParc,infoParc]=LoadVolumetricData([strcat(params.timepoint,'mo_Parc_Reg_Head')], ...
        fullfile('/Volumes/Extreme SSD/', strcat('mri/registered/UNC_to_NeuroDev/No Mask/', params.timepoint, 'mo')), ...
        'nii.gz');
end

% load mesh data
brainmeshLeft = load(strcat('/Volumes/Extreme SSD/mri/meshes/leftHemisphereMesh', params.timepoint,'mo.mat')); % loads nodes and faces
brainmeshRight = load(strcat('/Volumes/Extreme SSD/mri/meshes/rightHemisphereMesh', params.timepoint,'mo.mat'));

%% Search for task files 
matchingCortexFiles = testFuncs.getAgeTaskCortexDataFiles(params);
matchingParcelFiles = testFuncs.getAgeTaskParcelDataFiles(params);


%get count of files to be averaged
numParticipants = length(matchingCortexFiles);

groupHboSurfData = cell(params.nBlocks, 1);  % {1=baseline, 2–6=Blocks 1–5}
for i = 1:params.nBlocks
    groupHboSurfData{i} = {};  % inner cell array to store participant-level volumes
end

groupHbrSurfData = cell(params.nBlocks, 1);  % {1=baseline, 2–6=Blocks 1–5}
for i = 1:params.nBlocks
    groupHbrSurfData{i} = {};  % inner cell array to store participant-level volumes
end

expectedLength = []; % will hold length of surfData(:)
pendingNaNs = {};  % {blockIdx, iChrom}

for nsub = 1:length(matchingCortexFiles)

    [~, name, ~] = fileparts(matchingCortexFiles{nsub});
    fprintf('\nAnalysing file %d: %s\n', nsub, name);

    paramsFile = params;

    load(matchingCortexFiles{nsub});
    load(matchingParcelFiles{nsub});

    %find number of trials in each block
    nTrialsPerBlock = zeros(5, 1); 
    % Loop over 5 trial blocks: [1–5], [6–10], ..., [21–25]
    for iBlock = 1:5
        startTrial = (iBlock - 1) * 5 + 1;
        endTrial = iBlock * 5;
        
        % Logical mask of trials in this block
        inBlock = parcelData.trialNumbers >= startTrial & parcelData.trialNumbers <= endTrial;
        
        % Count how many are present
        nTrialsPerBlock(iBlock) = sum(inBlock);
    end
    blockIncluded = nTrialsPerBlock > 2;

    % load info for Jacobian matrix
    paramsFile.capName = analysisTools.getInfantCap(capCSV, capNames, paramsFile.timepoint, matchingCortexFiles{nsub}); %needs capName for jacobian name
    jacobFName = strcat('A_', paramsFile.capName, '_on_HD_Mesh_', paramsFile.timepoint, 'mo.mat'); %needs jacobian name to load it!
    infoJacob=load(fullfile(jacobianDir, jacobFName),'info'); %load info from the Jacobian

    for iChrom = 1:2
        chromData = cortexHbPeak{iChrom};  % cell array of block-level volumes
        nBlocksAvailable = numel(chromData);
        
        % Accumulate across available blocks (omit NaNs if desired)
        voxelSum = 0;
        voxelCount = 0;
        
        for blockIdx = 2:nBlocksAvailable
            included = blockIncluded(blockIdx-1);
            thisBlock = chromData{blockIdx};
            if isempty(thisBlock) || all(isnan(thisBlock(:)))

                if isempty(expectedLength)
                    % Defer assignment
                    pendingNaNs{end+1,1} = blockIdx;
                    pendingNaNs{end,2} = iChrom;
                else
                    placeholder = nan(expectedLength, 1);
                    if iChrom == 1
                        groupHboSurfData{blockIdx}{end+1} = placeholder;
                    else
                        groupHbrSurfData{blockIdx}{end+1} = placeholder;
                    end
                end
                continue

            end

            volData = Good_Vox2vol(thisBlock, infoJacob.info.tissue.dim);
            surfData = affine3d_img(volData, infoJacob.info.tissue.dim, infoParc);
            surfData(surfData == 0) = NaN; %remove 0s, to prevent biasing t-tests

            if isempty(expectedLength)
                expectedLength = numel(surfData(:));

                % fill any deferred NaNs
                for k = 1:size(pendingNaNs,1)
                    blk = pendingNaNs{k,1};
                    chrom = pendingNaNs{k,2};
                    placeholder = nan(expectedLength, 1);
                    if chrom == 1
                        groupHboSurfData{blk}{end+1} = placeholder;
                    else
                        groupHbrSurfData{blk}{end+1} = placeholder;
                    end
                end
            end
            
            if ~exist('surfD1', 'var')
                % save original dimensions of surface data for reshaping
                surfD1 = size(surfData, 1);
                surfD2 = size(surfData, 2);
                surfD3 = size(surfData, 3);
            end

            if iChrom == 1
                groupHboSurfData{blockIdx}{end+1} = surfData(:);
            else
                groupHbrSurfData{blockIdx}{end+1} = surfData(:);
            end


            if blockIdx == 2 && included == 1
                baseBlock = chromData{1};
                volData = Good_Vox2vol(baseBlock, infoJacob.info.tissue.dim);
                surfData = affine3d_img(volData, infoJacob.info.tissue.dim, infoParc);
                surfData(surfData == 0) = NaN; %remove 0s, to prevent biasing t-tests
                if iChrom == 1
                    groupHboSurfData{blockIdx-1}{end+1} = surfData(:);
                else
                    groupHbrSurfData{blockIdx-1}{end+1} = surfData(:);
                end

            end

        end

    end

end

%% HbO t-testing

% define pairs to be tested
% 1: baseline ---- 2: Hab1 ---- 3: Hab2 ---- 4: Hab3 ---- 5: Nov ---- 6:Post-test
testPairs = [2 1; 2 4; 2 5; 5 4]; 

hboTests = cell(length(testPairs), 1);

groupData = groupHboSurfData;

%perform test for each pair
for t = 1:size(testPairs, 1)

    fprintf(strcat('Testing condition pair ', num2str(t), '\n'))

    iA = testPairs(t, 1);
    iB = testPairs(t, 2);

    % ensure same number of subjects in both groups
    % should be fine but just in case
    nA = numel(groupData{iA});
    nB = numel(groupData{iB});
    nCommon = min(nA, nB);  % handle mismatch

    if nCommon < 3
        fprintf('Skipping test %d vs %d: not enough overlapping subjects.\n', iA, iB);
        continue
    end

    % create 2D data arrays - [voxels x subjects]
    dataA = cell2mat(groupData{iA}(1:nCommon));
    dataB = cell2mat(groupData{iB}(1:nCommon));

    % identify and exclude troublesome voxels across subjects
    validPairs = ~isnan(dataA) & ~isnan(dataB); 
    nValidPairs = sum(validPairs, 2);  
    validVoxels = nValidPairs >= nCommon./4; % include voxels with ≥ 25% valid pairs
    
    dataA_valid = dataA(validVoxels, :); 
    dataB_valid = dataB(validVoxels, :);
    [~, ~, ~, stats] = ttest(dataA_valid', dataB_valid'); %transpose so it's subjects x voxels, just for t-testing
    
    fullTstat = zeros(size(dataA, 1), 1);
    fullTstat(validVoxels) = stats.tstat;
    tMap = reshape(fullTstat, surfD1, surfD2, surfD3);

    % reshape and store t-map
    %tMap = reshape(stats.tstat, [], 1);
    %tMap(isinf(tMap) | isnan(tMap)) = 0; %put 0s back in for plotting
    %tMap = reshape(tMap, surfD1, surfD2, surfD3);

    hboTests{t} = tMap;

end




%% plot HbO surface data
tMaxVals = zeros(length(hboTests), 1);
for t = 1:length(hboTests)
    thisT = hboTests{t};
    tMaxVals(t) = max(abs(thisT(:)), [], 'omitnan');
end
paramsPlot.Scale = max(tMaxVals); % global max
paramsPlot.Scale = max(tMaxVals); % global max
paramsPlot.Th.P = (paramsPlot.Scale)./4; 
paramsPlot.Th.N = -paramsPlot.Th.P; 


for t = 1:length(hboTests)
    PlotInterpSurfMesh(hboTests{t}, brainmeshLeft, brainmeshRight, infoParc, paramsPlot);
end


%% HbR t-testing

% define pairs to be tested
% 1: baseline ---- 2: Hab1 ---- 3: Hab2 ---- 4: Hab3 ---- 5: Nov ---- 6:Post-test
testPairs = [2 1; 2 4; 2 5; 5 4]; 

hbrTests = cell(length(testPairs), 1);


groupData = groupHbrSurfData;

%perform test for each pair
for t = 1:size(testPairs, 1)

    fprintf(strcat('Testing condition pair ', num2str(t), '\n'))

    iA = testPairs(t, 1);
    iB = testPairs(t, 2);

    % ensure same number of subjects in both groups
    % should be fine but just in case
    nA = numel(groupData{iA});
    nB = numel(groupData{iB});
    nCommon = min(nA, nB);  % handle mismatch

    if nCommon < 3
        fprintf('Skipping test %d vs %d: not enough overlapping subjects.\n', iA, iB);
        continue
    end

    % create 2D data arrays - [voxels x subjects]
    dataA = cell2mat(groupData{iA}(1:nCommon));
    dataB = cell2mat(groupData{iB}(1:nCommon));

    % identify and exclude troublesome voxels across subjects
    validPairs = ~isnan(dataA) & ~isnan(dataB); 
    nValidPairs = sum(validPairs, 2);  
    validVoxels = nValidPairs >= nCommon./4; % include voxels with ≥ 25% valid pairs
    
    dataA_valid = dataA(validVoxels, :); 
    dataB_valid = dataB(validVoxels, :);
    [~, ~, ~, stats] = ttest(dataA_valid', dataB_valid'); %transpose so it's subjects x voxels, just for t-testing
    
    fullTstat = zeros(size(dataA, 1), 1);
    fullTstat(validVoxels) = stats.tstat;
    tMap = reshape(fullTstat, surfD1, surfD2, surfD3);

    hbrTests{t} = tMap;

end

%% plot HbR surface data
tMaxVals = zeros(length(hbrTests), 1);
for t = 1:length(hbrTests)
    thisT = hbrTests{t};
    tMaxVals(t) = max(abs(thisT(:)), [], 'omitnan');
end
paramsPlot.Scale = max(tMaxVals); % global max
paramsPlot.Scale = max(tMaxVals); % global max
paramsPlot.Th.P = (paramsPlot.Scale)./4; 
paramsPlot.Th.N = -paramsPlot.Th.P; 

for t = 1:length(hbrTests)
    PlotInterpSurfMesh(hbrTests{t}, brainmeshLeft, brainmeshRight, infoParc, paramsPlot);
end

%% HbT and HbD data construction
groupHbtSurfData = cell(params.nBlocks, 1);
groupHbdSurfData = cell(params.nBlocks, 1);

for iBlock = 1:params.nBlocks
    nSubs = min(numel(groupHboSurfData{iBlock}), numel(groupHbrSurfData{iBlock}));
    groupHbtSurfData{iBlock} = cell(1, nSubs);
    groupHbdSurfData{iBlock} = cell(1, nSubs);
    for s = 1:nSubs
        hbo = groupHboSurfData{iBlock}{s};
        hbr = groupHbrSurfData{iBlock}{s};
        if all(isnan(hbo)) || all(isnan(hbr))
            groupHbtSurfData{iBlock}{s} = nan(expectedLength, 1);
            groupHbdSurfData{iBlock}{s} = nan(expectedLength, 1);
        else
            groupHbtSurfData{iBlock}{s} = hbo + hbr;
            groupHbdSurfData{iBlock}{s} = hbo - hbr;
        end
    end
end

%% HbT t-testing
testPairs = [2 1; 2 4; 2 5; 5 4]; 
hbtTests = cell(length(testPairs), 1);
groupData = groupHbtSurfData;

for t = 1:size(testPairs, 1)

    fprintf(strcat('Testing condition pair ', num2str(t), '\n'))

    iA = testPairs(t, 1);
    iB = testPairs(t, 2);

    nA = numel(groupData{iA});
    nB = numel(groupData{iB});
    nCommon = min(nA, nB);

    if nCommon < 3
        fprintf('Skipping HbT test %d vs %d: not enough overlapping subjects.\n', iA, iB);
        continue
    end

    dataA = cell2mat(groupData{iA}(1:nCommon));
    dataB = cell2mat(groupData{iB}(1:nCommon));

    % identify and exclude troublesome voxels across subjects
    validPairs = ~isnan(dataA) & ~isnan(dataB); 
    nValidPairs = sum(validPairs, 2);  
    validVoxels = nValidPairs >= nCommon./4; % include voxels with ≥ 25% valid pairs
    
    dataA_valid = dataA(validVoxels, :); 
    dataB_valid = dataB(validVoxels, :);
    [~, ~, ~, stats] = ttest(dataA_valid', dataB_valid'); %transpose so it's subjects x voxels, just for t-testing
    
    fullTstat = zeros(size(dataA, 1), 1);
    fullTstat(validVoxels) = stats.tstat;
    tMap = reshape(fullTstat, surfD1, surfD2, surfD3);

    hbtTests{t} = tMap;
end

%% plot HbT surface data
tMaxVals = zeros(length(hbtTests), 1);
for t = 1:length(hbtTests)
    thisT = hbtTests{t};
    tMaxVals(t) = max(abs(thisT(:)), [], 'omitnan');
end
paramsPlot.Scale = max(tMaxVals); % global max
paramsPlot.Scale = max(tMaxVals); % global max
paramsPlot.Th.P = (paramsPlot.Scale)./4; 
paramsPlot.Th.N = -paramsPlot.Th.P; 

for t = 1:length(hbtTests)
    PlotInterpSurfMesh(hbtTests{t}, brainmeshLeft, brainmeshRight, infoParc, paramsPlot);
end

%% HbD t-testing
testPairs = [2 1; 2 4; 2 5; 5 4]; 
hbdTests = cell(length(testPairs), 1);
groupData = groupHbdSurfData;

for t = 1:size(testPairs, 1)

    fprintf(strcat('Testing condition pair ', num2str(t), '\n'))

    iA = testPairs(t, 1);
    iB = testPairs(t, 2);

    nA = numel(groupData{iA});
    nB = numel(groupData{iB});
    nCommon = min(nA, nB);

    if nCommon < 3
        fprintf('Skipping HbD test %d vs %d: not enough overlapping subjects.\n', iA, iB);
        continue
    end

    dataA = cell2mat(groupData{iA}(1:nCommon));
    dataB = cell2mat(groupData{iB}(1:nCommon));

    % identify and exclude troublesome voxels across subjects
    validPairs = ~isnan(dataA) & ~isnan(dataB); 
    nValidPairs = sum(validPairs, 2);  
    validVoxels = nValidPairs >= nCommon./4; % include voxels with ≥ 25% valid pairs
    
    dataA_valid = dataA(validVoxels, :); 
    dataB_valid = dataB(validVoxels, :);
    [~, ~, ~, stats] = ttest(dataA_valid', dataB_valid'); %transpose so it's subjects x voxels, just for t-testing
    
    fullTstat = zeros(size(dataA, 1), 1);
    fullTstat(validVoxels) = stats.tstat;
    tMap = reshape(fullTstat, surfD1, surfD2, surfD3);
    
    hbdTests{t} = tMap;
end

%% plot HbD surface data
tMaxVals = zeros(length(hbdTests), 1);
for t = 1:length(hbdTests)
    thisT = hbdTests{t};
    tMaxVals(t) = max(abs(thisT(:)), [], 'omitnan');
end
paramsPlot.Scale = max(tMaxVals); % global max
paramsPlot.Th.P = (paramsPlot.Scale)./4; 
paramsPlot.Th.N = -paramsPlot.Th.P; 

for t = 1:length(hbdTests)
    PlotInterpSurfMesh(hbdTests{t}, brainmeshLeft, brainmeshRight, infoParc, paramsPlot);
end