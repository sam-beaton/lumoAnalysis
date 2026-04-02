function checkSurfaceToDotAlignment(giiSurfFile, jacobInfo)
% checkSurfaceToDotAlignment
%
% INPUTS:
%   giiSurfFile : path to .surf.gii (HCP surface file)
%   jacobInfo   : jacob.info from your pipeline
%
% This function:
%   1) Loads surface vertices (MNI mm)
%   2) Builds voxel grid from DOT space
%   3) Compares bounding boxes + centres
%   4) Plots for visual sanity check
    
    %% --- Load surface ---
    g = gifti(giiSurfFile);
    vertices = double(g.vertices); % Nx3
    fprintf('\n--- SURFACE INFO ---\n');
    fprintf('Num vertices: %d\n', size(vertices, 1));
    
    % Bounding box (surface)
    surfMin = min(vertices, [], 1);
    surfMax = max(vertices, [], 1);
    surfCenter = mean(vertices, 1);
    fprintf('Surface bounding box (mm):\n');
    fprintf('  X: %.1f to %.1f\n', surfMin(1), surfMax(1));
    fprintf('  Y: %.1f to %.1f\n', surfMin(2), surfMax(2));
    fprintf('  Z: %.1f to %.1f\n', surfMin(3), surfMax(3));
    fprintf('Surface centre (mm): [%.1f %.1f %.1f]\n', surfCenter);
    
    %% --- Build DOT voxel grid ---
    xCoords = jacobInfo.tissue.dim.xv;
    yCoords = jacobInfo.tissue.dim.yv;
    zCoords = jacobInfo.tissue.dim.zv;
    
    nX = jacobInfo.tissue.dim.nVx;
    nY = jacobInfo.tissue.dim.nVy;
    nZ = jacobInfo.tissue.dim.nVz;
    
    [xGrid, yGrid, zGrid] = ndgrid(xCoords, yCoords, zCoords);
    voxelCoords = [xGrid(:), yGrid(:), zGrid(:)];
    
    % Bounding box (volume)
    volMin = min(voxelCoords, [], 1);
    volMax = max(voxelCoords, [], 1);
    volCenter = mean(voxelCoords, 1);
    fprintf('\n--- DOT VOLUME INFO ---\n');
    fprintf('Dims: [%d %d %d]\n', nX, nY, nZ);
    fprintf('Volume bounding box (mm):\n');
    fprintf('  X: %.1f to %.1f\n', volMin(1), volMax(1));
    fprintf('  Y: %.1f to %.1f\n', volMin(2), volMax(2));
    fprintf('  Z: %.1f to %.1f\n', volMin(3), volMax(3));
    fprintf('Volume centre (mm): [%.1f %.1f %.1f]\n', volCenter);
    
    %% --- Compare ---
    fprintf('\n--- COMPARISON ---\n');
    sizeDiff = (surfMax - surfMin) - (volMax - volMin);
    centerDiff = surfCenter - volCenter;
    fprintf('Size difference (surf - vol) [mm]: [%.1f %.1f %.1f]\n', sizeDiff);
    fprintf('Centre difference (surf - vol) [mm]: [%.1f %.1f %.1f]\n', centerDiff);
    
    %% --- Quick sanity check: sample vertices inside volume ---
    % Check how many vertices fall within volume bounds
    insideMask = vertices(:, 1) >= volMin(1) & vertices(:, 1) <= volMax(1) & ...
                 vertices(:, 2) >= volMin(2) & vertices(:, 2) <= volMax(2) & ...
                 vertices(:, 3) >= volMin(3) & vertices(:, 3) <= volMax(3);
    fprintf('Vertices inside volume: %.2f %%\n', 100 * mean(insideMask));
    
    %% --- Plot (subsample for speed) ---
    figure; hold on;
    
    % Subsample vertices
    nPlot = min(5000, size(vertices, 1));
    surfIdx = randperm(size(vertices, 1), nPlot);
    scatter3(vertices(surfIdx, 1), vertices(surfIdx, 2), vertices(surfIdx, 3), 5, 'b', 'filled');
    
    % Subsample voxel grid
    nVoxelPlot = min(5000, size(voxelCoords, 1));
    voxelIdx = randperm(size(voxelCoords, 1), nVoxelPlot);
    scatter3(voxelCoords(voxelIdx, 1), voxelCoords(voxelIdx, 2), voxelCoords(voxelIdx, 3), 5, 'r');
    
    legend({'Surface', 'DOT voxels'});
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Surface vs DOT space alignment');
    axis equal;
    grid on;
    view(3);

end