%% Script for Generating a pad file

%% NOTE
% When running this script, please DO NOT HIT THE GREEN TRIANGLE "RUN" BUTTON
% Please run individual code sections by highlighting and evalutating the highlighted selection
% Or, you can hit "ctrl + enter" on each section of code once the section's background turns beige
%
% SLB adapted from 

%% Pathing and output directory
clear all; close all; clc;

% Add toolboxes and relevant directories
addpath(genpath('/Users/sambe/Documents/MATLAB/toolboxes/NeuroDOT')); %neurodot toolbox

%% Define dir/file names here so all in one place and easier to change
% Set output directory
outputdir = ('/Volumes/G-DRIVE ArmorATD/imageRecon/neurodot/PADs'); %path to directory goes in between quotes
%Set cap name 
capName = 'GA00438'; % Create this yourself
%Set data file name - needs to match SD file for cap name
dataName = '/Volumes/G-DRIVE ArmorATD/dot/nirs/sub-062d/ses-06/nirs/sub-062d_ses-06_task-hand_run-01.nirs';
%% Change working dir to output path
if ~isfolder(outputdir)
    mkdir(outputdir);
end
cd(outputdir)

%% Load data and get optode locations

% load data that contains optode locations and wavelengths
Data = load(dataName, '-mat');

% create grid structure and place optode locations inside
    % grid structure should contain the following fields at a minimum
    % 3D source positions: grid.spos3
    % 3D detector positions: grid.dpos3
grid = struct;
grid.spos3 = Data.SD3D.SrcPos;
grid.dpos3 = Data.SD3D.DetPos;

% If available, place 2D optode positions in grid structure 
grid.spos2 = Data.SD.SrcPos(:, 1:2);
grid.dpos2 = Data.SD.DetPos(:, 1:2);


%% Create info structure (this is the pad file)

params.lambda = Data.SD3D.Lambda; % Ex: [750, 850]
params.mod = 'CW'; % Modulation type or frequency
params.CapName = capName; % set above
info = Generate_pad_from_grid(grid,params);


%% Visualize layout of pad file

% 2D layout - PlotCap
PlotCap(info)

% 3D layout - PlotCap
params_cap.dimension = '3D';
PlotCap(info, params_cap);view([-40,30])

%sbDotPlotCapWithLines(info, params_cap);view([-40,30])

% 3D layout - Draw_Foci
tpos=cat(1,info.optodes.spos3,info.optodes.dpos3); %SD positions
Ns=size(info.optodes.spos3,1); %Number of sources
Nd=size(info.optodes.dpos3,1); %Number of detectors
paramsFoci.color=cat(1,repmat([1,0,0],Ns,1),repmat([0,0,1],Nd,1)); % colors for sources (red) and detectors (blue)
paramsFoci.color(1,:) = [1 0.4 0.6]; %Pink for s1
paramsFoci.color(Ns+1,:) = [0.3010, 0.7450, 0.9330]; %Light blue for d1
figure;Draw_Foci_191203(tpos,paramsFoci);view([-40,30])


%% Make sanity check plots

% Visualize r2d by itself (hist)
figure;histogram(info.pairs.r2d,1000);xlabel('R_S_D');...
    ylabel('N_m_e_a_s');title('2D SD Separations');xlim([0 60])

% Visualize r3d by itself (hist)
figure;histogram(info.pairs.r3d,1000);xlabel('R_S_D');...
    ylabel('N_m_e_a_s');title('3D SD Separations');xlim([0 60])

% Visualize r3d by NN (scatter)
figure;scatter(info.pairs.r3d,info.pairs.NN);xlabel('R3d');ylabel('NN');
title([info.optodes.CapName, ', 3D SD sep by NN'], 'interpreter', 'none');

%% Match pad to data Section 
% Some data uses a different bookeeping strategy than NeuroDOT
% As such, it is important to make sure that whatever we make using
% NeuroDOT reflects the underlying data that we want to analyze
% This section will re-organize the pad file measurement list found in
% info.pairs to match the measurement list seen in the data
%%% Get Measurement lists from Pad and from Data

% Create temp pairs structure, this will get modified to match data
temp = info.pairs;

% Make measurement list from pad
pad_measList_before = [temp.Src, temp.Det, temp.WL];

% Make measurement list from data
data_measList = [Data.SD3D.MeasList(:,1), Data.SD3D.MeasList(:,2),...
    Data.SD3D.MeasList(:,3)];
data_measList = Data.SD3D.MeasList;

% Are measurement lists in the same order?
order_before = isequal(data_measList(:,[1,2,4]), pad_measList_before)

%%% Crop to match data
% Get order of measmurements in data
[~,Idx]=ismember(data_measList(:, [1,2,4]),pad_measList_before,'rows');

%Idx((length(Idx)/2)+1:end) = Idx((length(Idx)/2)+1:end) + length(Idx)/2;

% Re-order pairs structure
temp.r3d=[temp.r3d(Idx)];
temp.r2d=[temp.r2d(Idx)];
% temp.r2d = temp.r3d; %for sparse pads, uncomment this line - set r2d = to r3d
temp.NN = [temp.NN(Idx)];
temp.Src = [temp.Src(Idx)];
temp.Det = [temp.Det(Idx)];
temp.WL = [temp.WL(Idx);]; %Get correct size
temp.lambda = [temp.lambda(Idx)]; %Get correct size
temp.Mod = [temp.Mod(Idx,:),];

% Make sure measurements are identical
pad_measList_after = [temp.Src, temp.Det, temp.WL];
order_after = isequal(data_measList(:, [1,2,4]), pad_measList_after) %yes

% Sanity Check
% Visualize r3d by itself (hist)
figure;histogram(temp.r3d,1000);xlabel('R_S_D');ylabel('N_m_e_a_s');
title('SD Separations (Matched2Data)');xlim([0 60])

% 3D layout - PlotCap
params_cap.dimension = '3D';
PlotCap(info, params_cap);view([-40,30])

% 3D layout - Draw_Foci
tpos=cat(1,info.optodes.spos3,info.optodes.dpos3); %SD positions
Ns=size(info.optodes.spos3,1); %Number of sources
Nd=size(info.optodes.dpos3,1); %Number of detectors
paramsFoci.color=cat(1,repmat([1,0,0],Ns,1),repmat([0,0,1],Nd,1)); % colors for sources (red) and detectors (blue)
paramsFoci.color(1,:) = [1 0.4 0.6]; %Pink for s1
paramsFoci.color(Ns+1,:) = [0.3010, 0.7450, 0.9330]; %Light blue for d1
figure;Draw_Foci_191203(tpos,paramsFoci);view([-40,30])


%% Replace info.pairs with temp pairs and save pad file matched to data
info.pairs = temp;
save(['Pad_' info.optodes.CapName, '.mat'], 'info')