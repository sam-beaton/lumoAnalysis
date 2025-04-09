%% Setup and pathing
close all; clear all; 

% Add toolboxes and relevant directories
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

%% User-set file params
%%% Change where necessary
timePoint = '12'; %age of infant: '01', '06' or '12'
capName='GA00440'; %name of JSON file containing array info

%storage drive - easier than changing all names all the time
driveName = '/Volumes/G-DRIVE ArmorATD/';

jacobianDir=strcat(driveName, 'imageRecon/neurodot/Jacobians/');
meshDir = strcat(driveName, 'imageRecon/neurodot/Meshes/'); %folder containing meshes
outputDir=strcat(driveName, 'imageRecon/neurodot/workbench/'); %Output Directory for files

%%% Shouldn't need changing:
ldmeshname=[strcat('LD_Mesh_',timePoint,'mo')];     % LD Mesh name
hdmeshname=[strcat('HD_Mesh_',timePoint,'mo')];     % HD mesh name

%% Load necessary files
if ~isfolder(outputDir)
    mkdir(outputDir);
end
cd(outputDir);

% Load a Segmented Volume - YES
[maskSeg,infoSeg]=LoadVolumetricData([strcat(timePoint,'_0Months3T_head_segVol')],strcat(driveName, 'imageRecon/neurodot/Segmentations/', timePoint,'mo'),'nii');

% Load A matrix - ALREADY THERE
jacob = load([jacobianDir,'A_',capName,'_on_HD_Mesh_', timePoint, 'mo.mat'], 'info', 'A');

% Load LD mesh
%load([meshDir, ldmeshname,'.mat']);

% Load HD mesh
load([meshDir, hdmeshname,'.mat']);

% Load parcellation - YES
[maskParc,infoParc]=LoadVolumetricData([strcat(timePoint,'mo_Parc_Reg_Head')],strcat(driveName, 'mri/registered/UNC_to_NeuroDev/No Mask/', timePoint, 'mo'),'nii.gz');

% Extract 'dim' from PAD info as a separate variable for reg. to same space
% (not necessary but easier and consistent with conventions in Neurodot scripts)
dim = jacob.info.tissue.dim;

%% Visualise array
% params
pM.orientation='coord'; pM.Cmap.P='gray'; pM.EdgeColor='none'; pM.reg=0;
Ns=size(jacob.info.optodes.spos3,1);
Nd=size(jacob.info.optodes.dpos3,1);
paramsFoci.color=cat(1,repmat([1,0,0],Ns,1),repmat([0,0,1],Nd,1));
paramsFoci.color(1,:) = [1 0.4 0.6]; % pink for s
paramsFoci.color(Ns+1,:) = [0.3010, 0.7450, 0.9330]; %light blue for d
paramsFoci.radius = 2; %set radius of spheres to 2

% get relaxed optode poistions
tpos_relaxed = cat(1, jacob.info.optodes.spos3, jacob.info.optodes.dpos3);

% get mesh ready
m3=CutMesh(meshHD,intersect(find(meshHD.nodes(:,3)>0),find(meshHD.nodes(:,1)>0)));

%plot
PlotMeshSurface(m3,pM);Draw_Foci_191203(tpos_relaxed, paramsFoci);
view([-150,23]) %view mesh from behind with an angled top-down view

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

t1=affine3d_img(maskSeg,infoSeg,dim,eye(4)); % put anatomical volume in dim space

% Visualize Single SD pair (S1 D1)
keep=jacob.info.pairs.WL==2 & jacob.info.pairs.Src==1 & jacob.info.pairs.Det==1; % SD pair set here
foo=squeeze(jacob.A(keep,:));              % Single meas pair
fooV=Good_Vox2vol(foo',dim);
fooV=fooV./max(fooV(:));
fooV=log10(1e2.*fooV);                  % top 2 o.o.m.
pA.PD=1;pA.Scale=2;pA.Th.P=0;pA.Th.N=-pA.Th.P;
% PlotSlices(t1,dim,pA,fooV)

% Visualize Single SD pair (S20 D25)
keep=jacob.info.pairs.WL==2 & jacob.info.pairs.Src==20 & jacob.info.pairs.Det==25; % SD pair set here
foo=squeeze(jacob.A(keep,:));              % Single meas pair
fooV=Good_Vox2vol(foo',dim);
fooV=fooV./max(fooV(:));
fooV=log10(1e2.*fooV);                  % top 2 o.o.m.
pA.PD=1;pA.Scale=2;pA.Th.P=0;pA.Th.N=-pA.Th.P;
% PlotSlices(t1,dim,pA,fooV)

%% Visualise FFR

t1=affine3d_img(maskSeg,infoSeg,dim,eye(4)); % put anatomical volume in dim space

keep=(jacob.info.pairs.WL==2 & jacob.info.pairs.r3d<=60);
a=squeeze(jacob.A(keep,:));
iA=Tikhonov_invert_Amat(a,0.01,0.1);
iA=analysisTools.adaptedSmoothAmat(iA,dim,5); %5 = smoothing parameter
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

% sensitivity mask containing only GM
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

PlotSlices(t1,dim,paramParc,parcelsSens)


%% find and plot parts of parcellation not in GM mask 

t1_3 = t1.*(t1==3); %GM

% non-zero in parcel map but not GM map
mask = (parcelsSens ~= 0) & (t1_3 == 0);
% linear indices of those elements
idx = find(mask);
% Optionally, get the [row, col] subscripts if A and B are matrices
[row, col] = find(mask);

PlotSlices(t1,dim,pA,mask)



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
% Draw_Foci_191203(cat(1,jacob.info.optodes.spos3,jacob.info.optodes.dpos3),paramsFoci)  % Array
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