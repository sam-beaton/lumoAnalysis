% Edited to clear, close all previous work + add neurodot toolbox to path

clear; clc; close all

addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited function to save nii file in specified folder

%% Parameters for visualization
p.Cmap='jet'; p.Scale=5; p.Th.P=0; p.Th.N=-p.Th.P; p.PD=1; p.BG=[0,0,0];

timePoint = '12';
registeredImgsDir = strcat('/Users/sambe/mri/registered/UNC_to_NeuroDev/No Mask/', timePoint, 'mo');

%% Plot "_head" space together
[head, infoHead] = LoadVolumetricData( ...
    strcat('ANTS',timePoint,'-0Months3T_head'), ...
    '/Users/sambe/mri/NeurodevelopmentalMRI/Head', ...
    'nii.gz');
[brain, infoBrain] = LoadVolumetricData( ...
    strcat('ANTS',timePoint,'-0Months3T_head_brain'), ...
    '/Users/sambe/mri/NeurodevelopmentalMRI/Head', ...
    'nii.gz');

%PlotSlices(head, infoHead, [], brain) 

%% Load segmentations in '_head' space
% Either original (downloaded from NeuroDev MRI database or registered to
% same space

%Segmentation
[maskSeg, infoSeg] = LoadVolumetricData(strcat(timePoint,'mo_Seg_Reg_Head'), registeredImgsDir, 'nii.gz');
% PlotSlices(head, infoHead, [], maskSeg) 

%CSF
[maskCSF, infoCSF] = LoadVolumetricData(strcat(timePoint,'mo_CSF_Reg_Head'), registeredImgsDir, 'nii.gz');
% for i = size(maskCSF,1):-1:2
%     maskCSF(i, :, :) = maskCSF(i-1, :, :);
% end
% PlotSlices(head, infoHead, [], maskCSF) 

%GM
[maskGM, infoGM] = LoadVolumetricData(strcat(timePoint,'mo_GM_Reg_Head'), registeredImgsDir, 'nii.gz');
% PlotSlices(head, infoHead, [], maskGM) 

%WM
[maskWM, infoWM] = LoadVolumetricData(strcat(timePoint,'mo_WM_Reg_Head'), registeredImgsDir, 'nii.gz');
% PlotSlices(head, infoHead, [], maskWM) 

%Parcellation
[maskParc, infoParc] = LoadVolumetricData(strcat(timePoint,'mo_Parc_Reg_Head'), registeredImgsDir, 'nii.gz');
% PlotSlices(head, infoHead, [], maskParc) 

%% Binarize brain masks via logical indexing
mask_bin = zeros(size(head));
mask_bin(maskSeg == 10) = 1; %set csf to 1
mask_bin(maskSeg == 250) = 2; %set wm to 2
mask_bin(maskSeg == 150) = 3; %set gm to 3
% PlotSlices(mask_bin, infoHead);

%PlotSlices(head, infoHead, [], mask_bin) 

%% Get head mask

% Normalize
mask_norm = head./max(head(:));
% PlotSlices(mask_norm, infoHead)

% Extract
mask_head = ExtractHead(mask_norm);
% PlotSlices(mask_head,infoHead)

% Crop 
mask_crop = mask_head;
if strcmp(timePoint, '06') || strcmp(timePoint, '12')
    mask_crop(:,:,1:30) = 0; %at Z = 30
else
    mask_crop(:,:,1:25) = 0; %at Z = 35
end
% PlotSlices(mask_crop,infoHead)

%% Put Brain mask inside head mask
% initialize and set head to 4;
seg_mask = mask_crop;
seg_mask(seg_mask == 1) = 4;
% PlotSlices(seg_mask,infoHead)

% replace voxels in head mask with voxels from brain mask
seg_mask(mask_bin ~=0) = mask_bin(mask_bin ~=0);
% PlotSlices(seg_mask, infoHead); %visualize in grayscale

%% Save segmentation

outputDir = strcat('/Users/sambe/imageRecon/neurodot/Segmentations/', timePoint, 'mo');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

analysisTools.adaptedSaveVolumetricData(seg_mask, infoHead, strcat(timePoint,'_0Months3T_head_segVol'), outputDir, 'nii')
[mask,infoMask]=LoadVolumetricData(strcat(timePoint,'_0Months3T_head_segVol'), outputDir, 'nii');
PlotSlices(mask, infoMask, p)

