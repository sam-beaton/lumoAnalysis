%% Script to register a PAD file to a mesh generated from a segmented volume
% OVERVIEW
% 1. load mask and PAD
% 2. Split PAD into 3 pieces for registration
% 3. make a Low Density mesh of the mask
% 4. make a High Density mesh of the mask
% 5. use AlignMe to place the Pad on the LD mesh
% 6. run AlignMe on the LD-aligned Pad and the HD mesh
% 7. update the info file and name it specifically for this light model
% 8. prepare mesh for nirfast
% 9. set flags and run makeAnirfast
% 10. package and save A with info etc

%SLB 31/5/24
%Adapted from 'Basic_Light_Modeling_with_AlignMe_Tutorial.m'

% TO-DO:
% Automate entry of atlas fiducials Nz, Iz, M1, M2 - do youself for each mesh and use switch/case!
% Automate entry of S/Ds for each 'split' PAD - do youself for each array and use switch/case!
% Save the variables needed for visualisation (after constructing A) as don't want to have to do this each time!


%% Setup and pathing
close all; clear all

% Add toolboxes and relevant directories
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NIRFASTer')); %nirfast toolbox for meshing
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/Mesh2EEG')); %to generate fiducials for AlignMe
addpath(genpath('/Users/sambe/Documents/GitHubRepositories/nDotAnalysis')); %contains edited functions where necessary for use in image recon

% User-set file params
%%% Change where necessary
timePoint = '01'; %age of infant: '01', '06' or '12'
padname='GA00274'; %name of JSON file containing array info
arrayPositionAltered = 0;
arrayPosition = 'R';
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

% Load PAD file
load(strcat(padDir, '/Pad_', padname, '.mat')); 

%%Set parameters
%Parameter Initialization
p.Cmap='jet'; p.Scale=5; p.Th.P=0; p.Th.N=-p.Th.P; p.PD=1; p.BG=[0,0,0];
pM.orientation='coord'; pM.Cmap.P='gray'; pM.EdgeColor='none';
pS = pM; pS = rmfield(pS, 'orientation'); % make copy of pM, to be used with PlotSlices

%LD mesh params ---------------
param.facet_distance=5;    % Node position error tolerance at boundary
param.facet_size=3;        % boundary element size parameter
param.cell_size=5;         % Volume element size parameter
param.info=infoT1;
param.Offset=[0,0,0];
param.CheckMeshQuality=0;
param.Mode=1;              % make simple mesh with no region labels

%Array/pad file visualisation params -----------------------
tpos=cat(1,info.optodes.spos3,info.optodes.dpos3);
Ns=size(info.optodes.spos3,1);
Nd=size(info.optodes.dpos3,1);
paramsFoci.color=cat(1,repmat([1,0,0],Ns,1),repmat([0,0,1],Nd,1));
paramsFoci.color(1,:) = [1 0.4 0.6]; % pink for s
paramsFoci.color(Ns+1,:) = [0.3010, 0.7450, 0.9330]; %light blue for d
paramsFoci.radius = 2; %set radius of spheres to 2

%% Visualisation of loaded data
% Visualize the segmented mask
%PlotSlices(mask,infoT1,p)   

% Visualize PAD
%3D plot easier than 2D in this case to determine S/D lists
params_cap.dimension = '3D';
PlotCap(info, params_cap);view([-40,30]) %3D plot of pad
%PlotCap(info) %2D plot of pad

% Use 3D plot of pad to note where you want to split the pad up
switch padname
    case {'GA00274'}
        %sources - red
        S_1 = 1:18; %right
        S_2 = 19:36; %left
        %detectors - blue
        D_1 = 1:24; %right
        D_2  = 25:48; %left

    case {'GA00438_NF', 'GA00440_NF'}
        %sources - red
        S_1 = 10:27; %right
        S_2 = [1:9, 28:36]; %left
        %detectors - blue
        D_1 = 13:36; %right
        D_2 = [1:12, 37:48]; %left

    case {'GA00440', 'GA00438', 'GA00439'}
        %sources - red
        S_1 = 16:33; %right
        S_2 = [1:9, 37:45]; %left
        S_3 = [10:15, 34:36]; %frontal
        %detectors - blue
        D_1 = 21:44; %right
        D_2  = [1:12, 49:60]; %left
        D_3 = [13:20, 45:48]; %frontal

    case {'GA00370', 'GA00369', 'GA00351'}
        %sources - red
        S_1 = 13:27; %right
        S_2 = [1:9, 31:36]; %left
        S_3 = [10:12, 28:30]; %frontal
        %detectors - blue
        D_1 = 17:36; %right
        D_2  = [1:12, 41:48]; %left
        D_3 = [13:16, 37:40]; %frontal
end

% Split into separate pads

%%% Pad1, Right hand side ---------------------
params_cap.lambda = unique(info.pairs.lambda); % Ex: [750, 850]
params_cap.mod = 'CW'; % Modulation type or frequency
params_cap.CapName = 'Right_pad'; % Create this yourself
% Make pad1
tpos_Pad1 = cat(1,info.optodes.spos3(S_1,:),info.optodes.dpos3(D_1,:));
pad1 = info;
pad1.optodes.spos3 = pad1.optodes.spos3(S_1,:);
pad1.optodes.dpos3 = pad1.optodes.dpos3(D_1,:);
pad1.optodes.spos2 = pad1.optodes.spos2(S_1,:);
pad1.optodes.dpos2 = pad1.optodes.dpos2(D_1,:);
pad1 = Generate_pad_from_grid(pad1.optodes,params_cap);
pad1.pairs.r2d = pad1.pairs.r3d;
% Optode pos and visualization params for AlignMe
Ns_p1=size(pad1.optodes.spos3,1);
Nd_p1=size(pad1.optodes.dpos3,1);
paramsFoci_p1.color=cat(1,repmat([1,0,0],Ns_p1,1),repmat([0,0,1],Nd_p1,1));
paramsFoci_p1.color(1,:) = [1 0.4 0.6]; % pink for s1
paramsFoci_p1.color(Ns_p1+1,:) = [0.3010, 0.7450, 0.9330]; %light blue for d1
paramsFoci_p1.radius = 2; %set radius of spheres to 2

%%%Pad2, Left hand side -------------------------------
params_cap.CapName = 'Left_pad'; % Create this yourself
%Make pad2
tpos_Pad2 = cat(1,info.optodes.spos3(S_2,:),info.optodes.dpos3(D_2,:));
pad2 = info;
pad2.optodes.spos3 = pad2.optodes.spos3(S_2,:);
pad2.optodes.dpos3 = pad2.optodes.dpos3(D_2,:);
pad2.optodes.spos2 = pad2.optodes.spos2(S_2,:);
pad2.optodes.dpos2 = pad2.optodes.dpos2(D_2,:);
pad2 = Generate_pad_from_grid(pad2.optodes,params_cap);
pad2.pairs.r2d = pad2.pairs.r3d;
%optode pos and visualization params
Ns_p2=size(pad2.optodes.spos3,1);
Nd_p2=size(pad2.optodes.dpos3,1);
paramsFoci_p2.color=cat(1,repmat([1,0,0],Ns_p2,1),repmat([0,0,1],Nd_p2,1));
paramsFoci_p2.color(1,:) = [1 0.4 0.6]; % pink for s1
paramsFoci_p2.color(Ns_p2+1,:) = [0.3010, 0.7450, 0.9330]; %light blue for d1
paramsFoci_p2.radius = 2; %set radius of spheres to 2

%%%Pad3: Frontal -----------------------
if ~strcmp(gridname, 'GA00274') && ~contains(gridname, 'NF')
    params_cap.CapName = 'Front_pad'; % Create this yourself
    %Make pad3
    tpos_Pad3 = cat(1,info.optodes.spos3(S_3,:),info.optodes.dpos3(D_3,:));
    pad3 = info;
    pad3.optodes.spos3 = pad3.optodes.spos3(S_3,:);
    pad3.optodes.dpos3 = pad3.optodes.dpos3(D_3,:);
    pad3.optodes.spos2 = pad3.optodes.spos2(S_3,:);
    pad3.optodes.dpos2 = pad3.optodes.dpos2(D_3,:);
    pad3 = Generate_pad_from_grid(pad3.optodes,params_cap);
    pad3.pairs.r2d = pad3.pairs.r3d;
    %optode pos and visualization params
    Ns_p3=size(pad3.optodes.spos3,1);
    Nd_p3=size(pad3.optodes.dpos3,1);
    paramsFoci_p3.color=cat(1,repmat([1,0,0],Ns_p3,1),repmat([0,0,1],Nd_p3,1));
    paramsFoci_p3.color(1,:) = [1 0.4 0.6]; % pink for S
    paramsFoci_p3.color(Ns_p3+1,:) = [0.3010, 0.7450, 0.9330]; %light blue for D
    paramsFoci_p3.radius = 2; %set radius of spheres to 2
end

% Generate or load low density head mesh
if exist(strcat(meshDir, timePoint, 'mo/', ldmeshname, '.mat')) ~= 2
    %Generate LD mesh if one doesn't exist
    tic;meshLD=NirfastMesh_Region(mask,ldmeshname,param);toc
    % % Put coordinates back in true space
    meshLD.nodes=change_space_coords(meshLD.nodes,infoT1,'coord');
    PlotMeshSurface(meshLD,pM) %visualize in coordinate space
    view([135,30])  %this will give you a good view of where the center of the mesh is
    %save mesh for future use
    save(strcat(meshDir, ldmeshname, '.mat'), 'meshLD');
else
    %load LD mesh if it does exist
    load(strcat(meshDir, timePoint, 'mo/', ldmeshname, '.mat'));
    fprintf("Existing LD Mesh loaded\n");
end

% Visualize All 3 split pads in relation to LD mesh
% to make sure the whole array has been split correctly
% Right
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tpos_Pad1, paramsFoci_p1); %3D
view([90,0])

% Left
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tpos_Pad2, paramsFoci_p2); %3D
view([270,0])

% Front
if ~strcmp(gridname, 'GA00274') && ~contains(gridname, 'NF')
    PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tpos_Pad3, paramsFoci_p3); %3D
    view([270,0])
end

% Generate High Density Head Mesh or load existing mesh
close all;

if exist(strcat(meshDir, hdmeshname, '.mat')) ~=2
    %Create HD mesh if none exists in specified directory

    %If you get an error when running NirfastMesh_Region
    %try changing the mesh name and clearing your output directory of all files
    param.facet_distance=2.0;   % Node position error tolerance at boundary
    param.facet_size=0.8;       % boundary element size parameter
    param.cell_size=1.5;        % Volume element size parameter
    param.Mode=0;               % Set Mode=0 to make nirfast compliant mesh
    param.r0=5;                 % nodes outside of mask must be set to scalp==5;
    param.info=infoT1;          % Make sure info is infoT1 like with LD mesh
    param.Offset=[0,0,0];
    tic;meshHD=NirfastMesh_Region(mask,hdmeshname,param);toc
    
    %make copy of mesh that only contains nodes and elements for visualization purposes
    visMeshHD.nodes = meshHD.nodes;
    visMeshHD.elements = meshHD.elements;
    
    % Put coordinates back in coordinate space
    meshHD.nodes=change_space_coords(meshHD.nodes,infoT1,'coord'); %for actual mesh
    visMeshHD.nodes=change_space_coords(visMeshHD.nodes,infoT1,'coord'); %for mesh that's visualized

    %save mesh for future use
    save(strcat(meshDir, hdmeshname, '.mat'), 'meshHD');
    
    fprintf("HD Mesh generated and saved\n");
else
    %if HD mesh file exists, load it
    load(strcat(meshDir, hdmeshname, '.mat'));
    fprintf("HD Mesh loaded\n");
    %make copy of mesh that only contains nodes and elements for visualization purposes
    visMeshHD.nodes = meshHD.nodes;
    visMeshHD.elements = meshHD.elements;
    visMeshHD.nodes=change_space_coords(visMeshHD.nodes,infoT1,'coord'); %for mesh that's visualized
end 

% Find landmark coordinates
close all; %get rid of any stray plots hanging around
% Plot HD mesh; use cursor to select Nz, Iz etc and tner coordinates
% manually below

%PlotMeshSurface(visMeshHD,pM);view([70,60]) %Visualize in coordinate space

% Nz = Nasion position (indentation at the top of the nose approximately 
%   between the eyebrows
% Iz = Inion position (indentation at the back of the head approximately 
%   where the neck begins)
% M1: Left (as seen from above with nose on top) preauricular position 
%   (indentation in front of the top of the ear canal, dent between the 
%   upper edge of the tragus and the daith)
% M2: Right (as seen from above with nose on top) preauricular position 
%   (indentation in front of the top of the ear cannal, dent between the 
%    upper edge of the targus and the daith)
switch timePoint
    case '01'
        Nz = [0.4347 46.487 -3.9611]; 
        Iz = [0.3199 -42.3389 -62.237]; 
        M1 = [-38.248 -3.8748 -29.3093];
        M2 = [40.2975 -3.7481 -31.3124];
    case '06'
        Nz = [0.4637 51.661 -29.7167]; 
        Iz = [-0.6916 -72.4743 -63.1609]; 
        M1 = [-50.578 -20.8824 -37.2508];
        M2 = [51.3372 -18.2433 -43.1081]; 
    case '12'
        Nz = [-1.4448 59.805 -30.6696]; 
        Iz = [0.595 -63.38 -78.1934]; 
        M1 = [-54.377 -18.1152 -46.116];
        M2 = [56.3788 -14.045 -46.1];
end

% Find EEG (10-5) coordinates as inputs for AlignMe using Mesh2EEG/ComputeEEGPos
% Compile fiducial position matrix for input into function
m2eFiducs = [Nz;Iz;M1;M2];

% Use manually entered m2eFiducs to generate atlasFiducials
[EEGPts,EEGLab] = ComputeEEGPos(m2eFiducs, meshHD.elements, meshHD.nodes,1,0);

% Location of atlas fiducials                    Mesh2EEG output
atlasFiducials = [EEGPts(1,:); ... % Nasion --> EEGPts(1,:)
                  EEGPts(329,:); ... % Inion  EEGPts(329,:)
                  EEGPts(155,:); ... % LPA    EEGPts(155,:)
                  EEGPts(175,:); ... % RPA    EEGPts(175,:)
                  EEGPts(165,:)];    % Cz     EEGPts(165,:)

%% AlignMe Pad3 (frontal) (LD)
% this one first as easiest to align centrally wrt. nasion - can then use
% dy, dz, roty and rotz variables for other two pad arrays

if ~strcmp(gridname, 'GA00274') && ~contains(gridname, 'NF')
    % Create an instance of our custom DataStorage HANDLE class to store variables
    ds = DataStorage(); 
    
    % Input structure
    ds.dI.tpos = tpos_Pad3;   
    ds.dI.mesh = meshLD;
    ds.dI.pad = pad3; 
    ds.dI.pM = pM;
    ds.dI.paramsFoci = paramsFoci_p3; 
    ds.dI.Ns = Ns_p3;
    ds.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
    ds.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
    ds.dI.atlasFiducials = atlasFiducials;
    ds.dI.dy = -80;
    
    % Create an instance of your App Designer application, passing variable 'ds' into the app.  
    % The code is stuck at this line until your app closes, which destroys 'myapp'
    % But the data is assigned to ds variable in workspace
    % Run AlignMe
    myapp=AlignMe_2020b(ds);
    while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script
    
    % get relaxed optode positions from AlignMe
    tposNew = ds.dO.tpos2_relaxed; 
    tposNew_pad3 = tposNew;
    
    % visualize mesh with relaxed optodes
    pM.reg=0;
    PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew,paramsFoci_p3)
    view(270,25)
end

%% AlignMe Section, Pad1 (LD Mesh)  (RIGHT) 
% Create an instance of our custom DataStorage HANDLE class to store variables
ds = DataStorage(); 

% Input structure
ds.dI.tpos = tpos_Pad1;  
ds.dI.mesh = meshLD;
ds.dI.pad = pad1;          
ds.dI.pM = pM;
ds.dI.paramsFoci = paramsFoci_p1; 
ds.dI.Ns = Ns_p1;   
ds.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
ds.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
ds.dI.atlasFiducials = atlasFiducials;

% Create an instance of your App Designer application, passing variable 'ds' into the app.  
% The code is stuck at this line until your app closes, which destroys 'myapp'
% But the data is assigned to ds variable in workspace
% Run AlignMe
myapp=AlignMe_2020b(ds);
while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script

% get relaxed optode positions from AlignMe
tposNew = ds.dO.tpos2_relaxed;
tposNew_pad1 = tposNew;

% visualize mesh with relaxed optodes
pM.reg=0;
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew,paramsFoci_p1)
view(90,25)

%% AlignMe Pad2 (LD Mesh)   (LEFT)
% Create an instance of our custom DataStorage HANDLE class to store variables
ds = DataStorage(); 

% Input structure
ds.dI.tpos = tpos_Pad2;   
ds.dI.mesh = meshLD;
ds.dI.pad = pad2; 
ds.dI.pM = pM;
ds.dI.paramsFoci = paramsFoci_p2; 
ds.dI.Ns = Ns_p2;
ds.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
ds.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
ds.dI.atlasFiducials = atlasFiducials;

% Create an instance of your App Designer application, passing variable 'ds' into the app.  
% The code is stuck at this line until your app closes, which destroys 'myapp'
% But the data is assigned to ds variable in workspace
% Run AlignMe
myapp=AlignMe_2020b(ds);
while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script

% get relaxed optode positions from AlignMe
tposNew = ds.dO.tpos2_relaxed;
tposNew_pad2 = tposNew;

% visualize mesh with relaxed optodes
pM.reg=0;
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew,paramsFoci_p2)
view(270,25)

%% Relax optodes onto HD Mesh, pad3 (frontal)
% Again, relax this part first - should already be aligned, but for
% consistency with order of previous alignment on LD mesh
close all; %close any stray plots

if ~strcmp(gridname, 'GA00274') && ~contains(gridname, 'NF')

    % Create an instance of  DataStorage class
    ds_HD = DataStorage();
    
    % Input structure
    ds_HD.dI.tpos = tposNew_pad3;      
    ds_HD.dI.mesh = meshHD;
    ds_HD.dI.pad = pad3;            
    ds_HD.dI.pM = pM;
    ds_HD.dI.paramsFoci = paramsFoci_p3;
    ds_HD.dI.Ns = Ns_p3;
    ds_HD.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
    ds_HD.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
    ds_HD.dI.atlasFiducials = atlasFiducials;
    
    % Create an instance of your App Designer application,
    % passing variable 'ds' into the app.  
    % The code is stuck at this line until your app closes, which destroys 'myapp'
    % But the data is assigned to ds variable in workspace
    % Run AlignMe
    myapp = AlignMe_2020b(ds_HD);
    while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script
    
    % get relaxed optode positions from AlignMe
    tposNew_HD = ds_HD.dO.tpos2_relaxed;
    tposNew_HD_pad3 = tposNew_HD;
    
    % visualize mesh with relaxed optodes
    pM.reg=0;
    PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew_HD,paramsFoci_p3)
    view(90,25)
end

%% Relax optodes onto HD Mesh, pad1, right hand side
% Create an instance of  DataStorage class
ds_HD = DataStorage();

% Input structure
ds_HD.dI.tpos = tposNew_pad1;      
ds_HD.dI.mesh = meshHD;
ds_HD.dI.pad = pad1;            
ds_HD.dI.pM = pM;
ds_HD.dI.paramsFoci = paramsFoci_p1;
ds_HD.dI.Ns = Ns_p1;
ds_HD.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
ds_HD.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
ds_HD.dI.atlasFiducials = atlasFiducials;

% Create an instance of your App Designer application,
% passing variable 'ds' into the app.  
% The code is stuck at this line until your app closes, which destroys 'myapp'
% But the data is assigned to ds variable in workspace
% Run AlignMe
myapp = AlignMe_2020b(ds_HD);
while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script

% get relaxed optode positions from AlignMe
tposNew_HD = ds_HD.dO.tpos2_relaxed;
tposNew_HD_pad1 = tposNew_HD;

% visualize mesh with relaxed optodes
pM.reg=0;
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew_HD,paramsFoci_p1)
view(90,25)


%% Relax optodes onto HD Mesh, pad2, left hand side
% Create an instance of  DataStorage class
ds_HD = DataStorage();

% Input structure
ds_HD.dI.tpos = tposNew_pad2;   
ds_HD.dI.mesh = meshHD;
ds_HD.dI.pad = pad2;            
ds_HD.dI.pM = pM;
ds_HD.dI.paramsFoci = paramsFoci_p2;
ds_HD.dI.Ns = Ns_p2;
ds_HD.dI.photoPath = ''; %no participant photos for this tutorial, so leave as empty string
ds_HD.dI.imageType = '.jpeg'; %type of image file used, default is .jpeg
ds_HD.dI.atlasFiducials = atlasFiducials;

% Create an instance of your App Designer application,
% passing variable 'ds' into the app.  
% The code is stuck at this line until your app closes, which destroys 'myapp'
% But the data is assigned to ds variable in workspace
% Run AlignMe
myapp = AlignMe_2020b(ds_HD);
while isvalid(myapp); pause(0.1); end % Wait for app to close before continuing script

% get relaxed optode positions from AlignMe
tposNew_HD = ds_HD.dO.tpos2_relaxed;
tposNew_HD_pad2 = tposNew_HD;

% visualize mesh with relaxed optodes
pM.reg=0;
PlotMeshSurface(meshLD,pM);Draw_Foci_191203(tposNew_HD,paramsFoci_p2)
view(270,25)


%% Put relaxed pads back together, save PAD and check optode locations

%Find the numbers of sources and detectors, in order, as well as the
%(split) pade file in which they are contained. Then save this order, and
%concatenate to provide a full list of sources and detectors.
%This could probably be rewritten to run faster...

if ~strcmp(gridname, 'GA00274') && ~contains(gridname, 'NF')
    %sources
    newSourceInd = zeros([Ns, 1]);
    joinedSourceInd = [S_1, S_2, S_3];
    for iSrc = 1:Ns
        newSourceInd(iSrc) = find(joinedSourceInd == iSrc); %find index of iSrc'th source
    end
    %join source positions
    source_relaxed = cat(1,tposNew_HD_pad1(1:Ns_p1,:), tposNew_HD_pad2(1:Ns_p2,:), tposNew_HD_pad3(1:Ns_p3,:));
    %reorder according to original indices
    source_relaxed = source_relaxed(newSourceInd, :);
    
    %detectors
    newDetInd = zeros([Nd, 1]);
    joinedDetInd = [D_1, D_2, D_3];
    for iDet = 1:Nd
        newDetInd(iDet) = find(joinedDetInd == iDet);
    end
    %join detector positions
    detector_relaxed = cat(1,tposNew_HD_pad1(1+Ns_p1:end,:), tposNew_HD_pad2(1+Ns_p2:end,:), tposNew_HD_pad3(1+Ns_p3:end,:));
    %reorder according to original indices
    detector_relaxed = detector_relaxed(newDetInd, :);
    
    check_source = isequal(Ns, size(source_relaxed,1));
    check_detector = isequal(Nd, size(detector_relaxed,1));
    if check_source == 0 || check_detector == 0
        error('Sources or Detectors are not correct, please make sure you split the pad correctly')
    end
    tpos_relaxed = cat(1, source_relaxed, detector_relaxed);
else
    %sources
    newSourceInd = zeros([Ns, 1]);
    joinedSourceInd = [S_1, S_2];
    for iSrc = 1:Ns
        newSourceInd(iSrc) = find(joinedSourceInd == iSrc); %find index of iSrc'th source
    end
    %join source positions
    source_relaxed = cat(1,tposNew_HD_pad1(1:Ns_p1,:), tposNew_HD_pad2(1:Ns_p2,:));
    %reorder according to original indices
    source_relaxed = source_relaxed(newSourceInd, :);
    
    %detectors
    newDetInd = zeros([Nd, 1]);
    joinedDetInd = [D_1, D_2];
    for iDet = 1:Nd
        newDetInd(iDet) = find(joinedDetInd == iDet);
    end
    %join detector positions
    detector_relaxed = cat(1,tposNew_HD_pad1(1+Ns_p1:end,:), tposNew_HD_pad2(1+Ns_p2:end,:));
    %reorder according to original indices
    detector_relaxed = detector_relaxed(newDetInd, :);
    
    check_source = isequal(Ns, size(source_relaxed,1));
    check_detector = isequal(Nd, size(detector_relaxed,1));
    if check_source == 0 || check_detector == 0
        error('Sources or Detectors are not correct, please make sure you split the pad correctly')
    end
    tpos_relaxed = cat(1, source_relaxed, detector_relaxed);
end

%visualize
PlotMeshSurface(visMeshHD,pM);Draw_Foci_191203(tpos_relaxed,paramsFoci)
view(0,90) %dorsal view


% Update Pad info Structure 
info.optodes.spos3=tpos_relaxed(1:Ns,:);
info.optodes.dpos3=tpos_relaxed((Ns+1):end,:);
m=0;
for d=1:Nd
    for s=1:Ns
        m=m+1;
        info.pairs.r3d(m)=norm(tpos_relaxed(s,:)-tpos_relaxed(d+Ns,:)); %calculate new S-D distance
        info.pairs.r3d(m+(Ns*Nd))=info.pairs.r3d(m); %duplicate for second chrom
    end
end

% histogram of number of measurements per SD vs radius between SD pairs
% figure;histogram(info.pairs.r3d,1000);xlabel('R_S_D');ylabel('N_m_e_a_s');

% save pad file
save(['Pad_',hdmeshname,'_',gridname, '.mat'],'info') % CHANGE/ADD HERE

% PREPARE! --> mesh with grid array in same file set for NIRFAST
mesh=meshHD;
mesh=PrepareMeshForNIRFAST(mesh,[hdmeshname,'_',gridname],tpos_relaxed); % CHANGE/ADD POSITION TO PADNAME TO ALTER 
PlotMeshSurface(mesh,pM);Draw_Foci_191203(tpos_relaxed, paramsFoci);
view(0,90)

% One last visualization check...
% A portion of the mesh will be removed for this specific visualization
% This is so you can evaluate whether the optodes are placed too deeply within the mesh
m3=CutMesh(mesh,intersect(find(mesh.nodes(:,3)>0),find(mesh.nodes(:,1)>0)));
[Ia,Ib]=ismember(m3.nodes,mesh.nodes,'rows');Ib(Ib==0)=[];
m3.region=mesh.region(Ib);
PlotMeshSurface(m3,pM);Draw_Foci_191203(tpos_relaxed, paramsFoci);
view([-150,23]) %view mesh from behind with an angled top-down view

fprintf("Ready to calculate sensitivity profile\n\n");

%% Calculate Sensitivity Profile and save A matrix
% Set flags
flags.tag=[gridname,'_on_',hdmeshname];
flags.gridname=gridname;
flags.meshname=hdmeshname;
flags.head='info';
flags.info=infoT1;                                   % Your T1 info file
flags.gthresh=1e-5;                                  % Voxelation threshold in G
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
flags.srcnum=Ns;
% Get affine matrix that can be used to transform FROM participant space TO MNI space
if isfield(ds.dO, 'affineTform') %if mesh scaled, affineTform field will exist, save it to workspace
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
