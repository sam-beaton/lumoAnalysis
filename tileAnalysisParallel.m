%% Initialise and load toolboxes
clear all, close all;
% Toolbox paths
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT-main')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing - remove? only needed if generating meshes
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/lumoAnalysis')); %contains edited functions where necessary for use in image recon

%% Pathing and parameters
% ---------- User-defined parameters ------------
timepoints = {'06'};%{'01', '06', '12'};
params.task = 'hand'; %'hand', 'fc1' or 'fc2'

%storage drive - easier than changing all names all the time
driveName = '/Volumes/Extreme SSD/';

%overarching directory containing .nirs files
params.parentDir = fullfile(driveName, 'dot');
%processing method (%suffix after 'preproc-' in derivatives folder)
params.preProcDir = '025LPF'; 
% directory for (statistical) outputs
params.outputDir = fullfile(params.parentDir, 'derivatives'); %Output Directory for files

% ---------- Constant parameters - no need to change -------------
% cap names corresponding to info in capCSV file
capNames = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/capNames.csv';
%cap info for each participant, in .csv format
capCSV = '/Users/sambe/Library/CloudStorage/OneDrive-King''sCollegeLondon/Documents/INDiGO_docs/cappingData.csv'; 
% maximum channel distance to analyse
params.maxChannelDistance = 40; % Frijia et al (2021)
% Block averaging
params.dtPre = 40; % start
params.dtAfter = 190; % finish

% ---------- Derivative parameters -------------
%data location task and age dependant
params.dataLoc = fullfile(params.parentDir, 'derivatives', strcat('preproc-', params.preProcDir));

%%
% Start a parallel pool (specify number of workers)
% parpool(length(timepoints));

% parfor iTime = 1:length(timepoints)
for iTime = 1:length(timepoints)

    % Create a local copy of params for each worker
    paramsLocal = params;
    paramsLocal.timepoint = timepoints{iTime};

    % Search for task files 
    matchingFiles = analysisTools.getAgeTaskNirsFiles(paramsLocal);
    
    % Run Analysis
    for nsub = 1:length(matchingFiles) %01m: 59; 06mo: ? ; 12mo: 25
    
        [~, name, ~] = fileparts(matchingFiles{nsub});
        fprintf('\nAnalysing file %d: %s\n', nsub, name);
    
        % ------ Reset file-specific parameters -----
        paramsFile = paramsLocal;
    
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
            
            % ---------- Get trial numbers ------------
            trialNumbers = analysisTools.getTrialNumbers(data.info);
            
            % ------- Calculate block averaged data ----------
            %[badata, ~, ~, ~, tKeep] = analysisTools.adaptedBlockAverage(lmdata, paramsFile, data.info);
            
            % ------- View block averaged data ---------
            %analysisTools.viewBlockAveraged(badata, paramsFile);
            
            % ---------- Obtain block data for each channel ----------
            [channelBlockData, sourceNumbers, detectorNumbers, paramsFile] = analysisTools.getChannelBlockData(data.dc, data.info, paramsFile);
            
            % ---------- Checks and housekeeping ------------   
            %check trial numbers make sense - continue to next file if not
            if isempty(trialNumbers)
                fprintf('synchtypes does not match expected pattern. Skipping...\n');
                continue; % Skip to next iteration
            end
    
            if isempty(channelBlockData)
                continue; % Skip to next iteration - warning already in getChannelBlockData
            end
            
            % check all blocks were used for averaging i.e. weren't too close to 
            % ends of recording
            if isfield(paramsFile, 'blockRemoved')
                trialNumbers(paramsFile.blockRemoved) = []; %removes unused trial number (1st or last)
            end
    
            % ---------- Get tile-averaged block data -------------
            [tileBlockData, tileNumbers] = analysisTools.getTileBlockData(channelBlockData, sourceNumbers, detectorNumbers, paramsFile);
    
            % ---------- Join channel data ____________
            channelData = struct;
            channelData.blockData = channelBlockData;
            channelData.trialNumbers = trialNumbers;
            channelData.sourceNumbers = sourceNumbers;
            channelData.detectorNumbers = detectorNumbers;
            channelData.capName = paramsFile.capName;
    
            % --------- Join tile data -----------
            tileData = struct;
            tileData.blockData = tileBlockData; %unfinished
            tileData.tileNumbers = tileNumbers; % unfinished
            tileData.trialNumbers = trialNumbers;
            tileData.capName = paramsFile.capName;
            
            % ---------- Save variables ----------
            analysisTools.saveFiles(paramsFile, matchingFiles{nsub}, [], [], [], channelData, tileData);
            
        catch
            fprintf(strcat('Could not run analysis - look into manually.\n'))
        end
    
        % ---------- Tidying up ----------
        % prevent stale variable carry-over:
        data = [];
        lmdata = [];
        channelBlockData = [];
        tileBlockData = [];
        tileNumbers = [];
        trialNumbers = [];
        close all; %incase plotting used
    end
end

% Shut down the parallel pool
% delete(gcp('nocreate'));

fprintf("\n=========================FINISHED PROCESSING===========================\n")

%% TESTING


