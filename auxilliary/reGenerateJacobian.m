% Script to re-generate the Jacobian, after first altering the parameters

% Setup and pathing
close all; clear all

% Add toolboxes and relevant directories
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/Mesh2EEG')); %to generate fiducials for AlignMe
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

% User-set file params
%%% Change where necessary
timePoint = '12'; %age of infant: '01', '06' or '12'
padname='GA00440'; %name of JSON file containing array info
arrayPositionAltered = 1;
arrayPosition = 'U';
driveName = '/Volumes/Extreme SSD/'; %storage drive - easier than changing all names all the time
meshDir = fullfile(driveName, 'imageRecon/neurodot/Meshes/');
outputDir = fullfile(driveName, 'imageRecon/neurodot/Jacobians/');
padDir = fullfile(driveName, 'imageRecon/neurodot/PADs/');
%meshDir = '/Volumes/G-DRIVE ArmorATD/imageRecon/neurodot/Meshes/'; %folder to save meshes to/load them from
%outputDir = '/Volumes/G-DRIVE ArmorATD/imageRecon/neurodot/Jacobians'; %Output Directory for files

%%% Shouldn't need changing:
ldmeshname=[strcat('LD_Mesh_',timePoint,'mo')];     % LD Mesh name
hdmeshname=[strcat('HD_Mesh_',timePoint,'mo')];     % HD mesh name

% Set name for *saving* pade and associated objects, based on position
if arrayPositionAltered == 1
    gridname = strcat(padname, '_', arrayPosition);
else
    gridname = padname;
end

if ~isfolder(meshDir)
    mkdir(meshDir);
end

if ~isfolder(outputDir)
    mkdir(outputDir);
end

% Load necessary files
cd(outputDir)

% Load a Segmented Volume
[mask,infoT1]=LoadVolumetricData([strcat(timePoint,'_0Months3T_head_segVol')],strcat(driveName, 'imageRecon/neurodot/Segmentations/', timePoint,'mo'),'nii');

% load pad on mesh
load(['Pad_',hdmeshname,'_',gridname, '.mat']);

%if HD mesh file exists, load it
load(strcat(meshDir, hdmeshname, '.mat'));
fprintf("HD Mesh loaded\n");

%prepare mesh for light modelling
mesh=meshHD;
tpos_relaxed = cat(1, info.optodes.spos3, info.optodes.dpos3);
mesh=PrepareMeshForNIRFAST(mesh,[hdmeshname,'_',gridname],tpos_relaxed);

%% Calculate Sensitivity Profile and save A matrix
% Set flags
flags.tag=[gridname,'_on_',hdmeshname];
flags.gridname=gridname;
flags.meshname=hdmeshname;
flags.head='info';
flags.info=infoT1;                                   % Your T1 info file
flags.gthresh=1e-3;                                  % Voxelation threshold in G
flags.voxmm=2;                                       % Voxelation resolution (mm)
flags.labels.r1='csf';                               % Regions for optical properties
flags.labels.r2='white';
flags.labels.r3='gray';
flags.labels.r4='bone';
flags.labels.r5='skin';
flags.op.lambda=[750,850];                           % Wavelengths (nm)
flags.op.mua_skin=[0.0170,0.0190];                   % Baseline absorption
flags.op.mua_bone=[0.0116,0.0139];
flags.op.mua_csf=[0.0040,0.0040];
flags.op.mua_gray=[0.0180,0.0192];
flags.op.mua_white=[0.0167,0.0208];
flags.op.musp_skin=[0.74,0.64];                      % Baseline reduced scattering coeff
flags.op.musp_bone=[0.94,0.84];
flags.op.musp_csf=[0.3,0.3];
flags.op.musp_gray=[0.8359,0.6726];
flags.op.musp_white=[1.1908,1.0107];
flags.op.n_skin=[1.4,1.4];                           % Index of refraction
flags.op.n_bone=[1.4,1.4];
flags.op.n_csf=[1.4,1.4];
flags.op.n_gray=[1.4,1.4];
flags.op.n_white=[1.4,1.4];
flags.srcnum=size(info.optodes.spos3,1);
% Get affine matrix that can be used to transform FROM participant space TO MNI space
if exist('ds', 'var') && isfield(ds.dO, 'affineTform') %if mesh scaled, affineTform field will exist, save it to workspace
    affine_Subj2MNI = [ds.dO.affineTform, zeros(3,1)];
    save('affine_matrix_Subject_to_MNI.mat', 'affine_Subj2MNI')
else %otherwise, set affine_Subj2MNI to eye(4)
    affine_Subj2MNI = eye(4);
end% Number of sources
flags.t4=affine_Subj2MNI;                            % Affine matrix for going from subject-specific space to MNI space
flags.t4_target='MNI';                               % string
flags.makeA=1;                                       % don't make A, just make G
flags.Hz=0;
if flags.Hz, flags.tag = [flags.tag,'FD']; end

% Run makeAnirfast to get sensitivity matrix
Ti=tic;[A,dim,Gsd]=makeAnirfaster(mesh,flags); % size(A)= [Nwl, Nmeas, Nvox]
disp(['<makeAnirfast took ',num2str(toc(Ti))])

% %Gsd vs Rsd provides a simulated light fall-off
% figure;semilogy(info.pairs.r3d(1:(Ns*Nd)),Gsd,'*');
% legend({num2str(flags.op.lambda')})
% xlabel('R_S_D [mm]');ylabel('G_S_D');xlim([0,100])


% Package data and save A
[Nwl,Nmeas,Nvox]=size(A); %Nwl=#wavelengths; Nmeas=##channels; Nvox=#good channels in voxel interpolation of light sensitivity mask
A=reshape(permute(A,[2,1,3]),Nwl*Nmeas,Nvox); %converts from 3D to 2D matrix by concatenating wrt. channels i.e. is now a 2D (Nwl*Nmeas)xNvox matrix

% Place spatial information about light model in info.tissue structure
info.tissue.dim=dim; %contains bookkeeping/metadata about sensitivity mask in voxel space
info.tissue.affine=flags.t4;
info.tissue.infoT1=infoT1;
info.tissue.affine_target='MNI';
info.tissue.flags=flags;

save(['A_',flags.tag,'.mat'],'A','info','-v7.3') %save A

fprintf("A matrix reshaped and saved\n");

