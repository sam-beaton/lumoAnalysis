%% Setup and pathing
close all; clear all; 

% Add toolboxes and relevant directories
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

%% User-set file params
%%% Change where necessary
timePoint = '06'; %age of infant: '01', '06' or '12'
capName='GA00370_MisLeft'; %name of JSON file containing array info
jacobianDir=strcat('/Volumes/G-DRIVE ArmorATD/imageRecon/neurodot/Jacobians/', capName, '_', timePoint, 'mo/');
meshDir = strcat('/Volumes/G-DRIVE ArmorATD/imageRecon/neurodot/Meshes/', timePoint, 'mo/'); %folder containing meshes
outputDir=strcat('/Volumes/G-DRIVE ArmorATD/imageRecon/neurodot/workbench/', capName, '_', timePoint, 'mo'); %Output Directory for files

%%% Shouldn't need changing:
ldmeshname=[strcat('LD_Mesh_',timePoint,'mo')];     % LD Mesh name
hdmeshname=[strcat('HD_Mesh_',timePoint,'mo')];     % HD mesh name

%% Load necessary files
if ~isfolder(outputDir)
    mkdir(outputDir);
end
cd(outputDir);

% Load a Segmented Volume
[maskSeg,infoSeg]=LoadVolumetricData([strcat(timePoint,'_0Months3T_head_segVol')],strcat('/Volumes/G-DRIVE ArmorATD/imageRecon/neurodot/Segmentations/', timePoint,'mo'),'nii');

% Load PAD file
load([jacobianDir, 'Pad_HD_Mesh_', timePoint, 'mo_', capName, '.mat']); %loads as 'info'

% Load A matrix
load([jacobianDir,'A_',capName,'_on_HD_Mesh_', timePoint, 'mo.mat']);

% Load LD mesh
load([meshDir, ldmeshname,'.mat']);

% Load HD mesh
load([meshDir, hdmeshname,'.mat']);

% Load parcellation
[maskParc,infoParc]=LoadVolumetricData([strcat(timePoint,'mo_Parc_Reg_Head')],strcat('/Users/sambe/mri/registered/UNC_to_NeuroDev/No Mask/', timePoint, 'mo'),'nii.gz');

% Extract 'dim' from PAD info as a separate variable for reg. to same space
% (not necessary but easier and consistent with conventions in Neurodot scripts)
dim = info.tissue.dim;

%% Set params for plotting etc
% PlotSlices params ---------
pA.PD=1;
pA.Scale=2;
pA.Th.P=0;
pA.Th.N=-pA.Th.P;

% HD head surface mesh params
pA2.fig_handle=gca;
pA2.FaceAlpha=0;
pA2.EdgeAlpha=0.5;
pA2.EdgeColor=[1,1,1].*0.25;

% cortical surface mesh visualisation params
pA0l=struct;
pA0l.fig_handle=gca;
pA0l.FaceColor=[1,1,1].*0.5;
pA0l.EdgeColor='none';
pA0l.AmbientStrength=0.25;
pA0l.DiffuseStrength=0.25;
pA0l.SpecularStrength=0.025;

% PlotMeshSurface params ------------
pM.orientation='coord'; pM.Cmap.P='gray'; pM.EdgeColor='none';
% PlotSlice params ------------ 
pS = pM; pS = rmfield(pS, 'orientation'); % make copy of pM, to be used with PlotSlices

%LD mesh params ---------------
paramLD.facet_distance=5;    % Node position error tolerance at boundary
paramLD.facet_size=3;        % boundary element size parameter
paramLD.cell_size=5;         % Volume element size parameter
paramLD.info=infoSeg;
paramLD.Offset=[0,0,0];
paramLD.CheckMeshQuality=0;
paramLD.Mode=1;              % make simple mesh with no region labels

% HD mesh params --------------
paramHD.facet_distance=2.0;   % Node position error tolerance at boundary
paramHD.facet_size=0.8;       % boundary element size parameter
paramHD.cell_size=1.5;        % Volume element size parameter
paramHD.info=infoSeg;          % Make sure info is infoT1 like with LD mesh
paramHD.Offset=[0,0,0];
paramHD.Mode=0;               % Set Mode=0 to make nirfast compliant mesh
paramHD.r0=5;                 % nodes outside of mask must be set to scalp==5;

close all;

%% Visualize aspects of sensitivity profile

% dim.center(1) = dim.center(1)+2; %shift over by 1 voxel (2mm)
info.tissue.dim = dim;

t1=affine3d_img(maskSeg,infoSeg,dim,eye(4)); % put anatomical volume in dim space

% Visualize Single SD pair (S1 D1)
keep=info.pairs.WL==2 & info.pairs.Src==1 & info.pairs.Det==1; % SD pair set here
foo=squeeze(A(keep,:));              % Single meas pair
fooV=Good_Vox2vol(foo',dim);
fooV=fooV./max(fooV(:));
fooV=log10(1e2.*fooV);                  % top 2 o.o.m.
pA.PD=1;pA.Scale=2;pA.Th.P=0;pA.Th.N=-pA.Th.P;
% PlotSlices(t1,dim,pA,fooV)

% Visualize Single SD pair (S20 D25)
keep=info.pairs.WL==2 & info.pairs.Src==20 & info.pairs.Det==25; % SD pair set here
foo=squeeze(A(keep,:));              % Single meas pair
fooV=Good_Vox2vol(foo',dim);
fooV=fooV./max(fooV(:));
fooV=log10(1e2.*fooV);                  % top 2 o.o.m.
pA.PD=1;pA.Scale=2;pA.Th.P=0;pA.Th.N=-pA.Th.P;
% PlotSlices(t1,dim,pA,fooV)

%% Visualise FFR
keep=(info.pairs.WL==2 & info.pairs.r3d<=60);
a=squeeze(A(keep,:));
iA=Tikhonov_invert_Amat(a,0.01,0.1);
iA=smooth_Amat(iA,dim,5); %5 = smoothing parameter
ffr=makeFlatFieldRecon(a,iA);

fooV=Good_Vox2vol(ffr,dim);
fooV=fooV./max(fooV(:));
pA.PD=1;pA.Scale=1;pA.Th.P=5e-2;pA.Th.N=-pA.Th.P;
PlotSlices(t1,dim,pA,fooV)


%% Find intersection of light model and GM
% find indices of voxels imaged by array, GM voxels, and common indices to
% both
sensMask = find(fooV > 0.05); %find indices of light coverage mask
indGM = find(t1 == 3); %find indices of GM in anatomical volume
gmLight = intersect(sensMask, indGM); %find common indices i.e. GM covered by array

% obtain GM space covered by array
gmIndBin  = zeros([size(t1)]);
gmIndBin(gmLight) = 1; %create binary mat. to multiply by FooV
GMfooV = gmIndBin.*fooV; %multiply to remove non-GM tissues from overlay

PlotSlices(t1,dim,pA,GMfooV)

%% Find intersection of light model and parcellation
% transform Parcellated GM volume to dim space
regMaskParc=affine3d_img(maskParc,infoSeg,dim,eye(4)); 

%find common indices i.e. parcellated GM covered by sensitivity mask
indRegMaskParc = find(regMaskParc ~= 0);
parcLight = intersect(sensMask, indRegMaskParc); 
 
% obtain GM space covered by array
indRegParcBin  = zeros([size(t1)]);
indRegParcBin(parcLight) = 1; %create binary mat. to multiply by FooV

% sensitivity mask containing only GM/parcellated voxels
fooVParcOnly = indRegParcBin.*fooV;
PlotSlices(t1,dim,pA,fooVParcOnly)

% parcels overlapping with sensitivity mask
parcFooVOnly = indRegParcBin.*regMaskParc;
pA.cbmode = 0;
PlotSlices(t1,dim,pA,parcFooVOnly)

%% Define parcels which mask is sensitive to
% ================================== FINE TO HERE! ===========================================
parcelsIncl = unique(regMaskParc); %find all parcel numbers in a vector
parcelsIncl(1) = []; %remove 0 as it corresponds to non-GM space
%initialise mask of parcels which array is sensitive to:
parcelsSens = zeros(size(regMaskParc)); 
for iParc = length(parcelsIncl):-1:1
    %check proportion of parcel voxels covered by array
    if length(find(parcFooVOnly == parcelsIncl(iParc))) / length(find(regMaskParc == parcelsIncl(iParc))) >= 0.5
        %if >= 50%, include parcel in mask
        parcelsSens(find(regMaskParc == parcelsIncl(iParc))) = parcelsIncl(iParc);
    else
        %otherwise remove it from list of included parcels
        parcelsIncl(iParc)=[];
    end
end

paramParc=struct;
paramParc.TC = 1; %use true color mapping
paramParc.cbmode = 0; %use custom colorbar limits
paramParc.Cmap.P = colorcube; %use custom colormap
paramParc.PD = 1; %values in image positive-definite

sbDotPlotParcelSlices(t1,dim,paramParc,parcelsSens)

% %%
% Smesh = vol2surf_mesh(meshParc, maskParc, dim);
% 
% %% Visualize alignment of LD mesh, HD mesh, array, and cortex
% % create GM mesh
% maskGM = zeros(size(maskSeg)); maskGM(maskSeg == 3) = 1;
% 
% paramGM = paramHD;
% paramGM.facet_size = 0.15; 
% paramGM.cell_size = 0.2;
% paramGM.facet_angle = 30;
% paramGM.mode = 1;
% tic;meshGM=NirfastMesh_Region(maskGM,strcat('HD_Mesh_GM_', timePoint,'mo'), paramGM);toc;
% meshGM.nodes=change_space_coords(meshGM.nodes,infoSeg,'coord');
% PlotMeshSurface(meshGM,pM) 
% visMeshGM.nodes = meshGM.nodes;
% visMeshGM.elements = meshGM.elements;
% PlotMeshSurface(visMeshGM,pM);view([70,60])
% 
% %  View alignment
% pA0l=struct;
% figure('Color','k','Position',[500,100,1050,1000])
% pA0l.fig_handle=gca;
% pA0l.FaceColor=[1,1,1].*0.5;pA0l.EdgeColor='none';
% pA0l.AmbientStrength=0.25;pA0l.DiffuseStrength=0.25;
% pA0l.SpecularStrength=0.025;
% PlotMeshSurface(meshGM,pA0l)                 % Cortical Surface
% pA2.fig_handle=gca;
% pA2.FaceAlpha=0;
% pA2.EdgeAlpha=0.5;
% pA2.EdgeColor=[1,1,1].*0.25;
% PlotMeshSurface(meshHD,pA2)                     % LD mesh
% set(gca,'Color','k')
% view([90,0]) % display in Lateral view
% Draw_Foci_191203(cat(1,info.optodes.spos3,info.optodes.dpos3),paramsFoci)  % Array
% axis off
% 
% %% visualise parcellation in mesh space
% paramParc = paramHD;
% paramParc.Mode = 1;
% paramParc.facet_angle = 30;
% paramParc.cell_size = 0.1;
% tic;meshParc=NirfastMesh_Region(maskParc,'HD_Mesh_Parc_12mo', paramParc);toc;
% visMeshParc.nodes = meshParc.nodes;
% visMeshParc.elements = meshParc.elements;
% pMP = pM;
% pMP.Cmap = 'hsv';
% PlotMeshSurface(meshParc,pMP);view([70,60])
% 
% [f,v] = vol2isomesh(parc, infoParc);