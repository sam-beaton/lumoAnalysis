function volParc = surfaceParcToVolume(giiSurfFile, labelFile, jacobInfo)
% surfaceParcToVolume
% Converts a surface parcellation to a DOT volumetric label image
%
% INPUTS:
%   giiSurfFile : path to .surf.gii (HCP surface, FSL voxel mm space)
%   labelFile   : path to .label.gii or .func.gii with parcel labels per vertex
%   jacobInfo   : jacob.info from your pipeline
%
% OUTPUT:
%   volParc     : 3D label volume in DOT voxel space, same dims as DOT grid

%% --- Load surface vertices and labels ---
g = gifti(giiSurfFile);
vertices = double(g.vertices);   % Nx3, in FSL mm space

L = gifti(labelFile);
labels = double(L.cdata);        % Nx1, integer parcel labels

%% --- Transform vertices from FSL mm → MNI mm ---
% jacob.info.tissue.affine is 4x4, maps from voxel index → MNI mm
% But surface is already in FSL mm, not voxel index.
% We need the T1 vox2world matrix to go: FSL mm → voxel → MNI
%
% Option A: if you have the T1 nii, extract sform (recommended)
%   nii = load_nii('T1.nii.gz');
%   vox2world = nii.hdr.hist.sform_mat; % 4x4
%
% Option B: use affine directly if it maps FSL mm → MNI mm
% Try this first and check with checkSurfaceToDotAlignment afterwards:
affine = jacobInfo.tissue.affine;  % 4x4

% Homogeneous coordinates
verticesH = [vertices, ones(size(vertices,1), 1)]';  % 4xN
verticesMni = (affine * verticesH)';                  % Nx4
verticesMni = verticesMni(:, 1:3);                    % Nx3, now in MNI mm

fprintf('Transformed vertices centre: [%.1f %.1f %.1f]\n', mean(verticesMni));

%% --- Build DOT voxel coordinate vectors ---
xv = jacobInfo.tissue.dim.xv(:);   % column vectors
yv = jacobInfo.tissue.dim.yv(:);
zv = jacobInfo.tissue.dim.zv(:);

nX = jacobInfo.tissue.dim.nVx;
nY = jacobInfo.tissue.dim.nVy;
nZ = jacobInfo.tissue.dim.nVz;

%% --- Map each vertex to its nearest voxel ---
% Use Gaussian-weighted parcel voting per voxel
sigma = 3;           % mm - Gaussian kernel width
radius = 6;          % mm - search radius (2*sigma)

% Initialise accumulator: nVoxels x nParcels (sparse would be better for large grids)
nParcels = max(labels);
voxWeights = zeros(nX, nY, nZ, nParcels, 'single');

fprintf('Mapping %d vertices to voxel grid...\n', size(verticesMni,1));

% Precompute: find bounding voxel indices within radius for each vertex
% Vectorised over vertices in chunks for memory efficiency
chunkSize = 5000;
nVerts = size(verticesMni, 1);

for vStart = 1:chunkSize:nVerts

    vEnd = min(vStart + chunkSize - 1, nVerts);
    vChunk = verticesMni(vStart:vEnd, :);   % chunk x 3
    lChunk = labels(vStart:vEnd);           % chunk x 1

    for v = 1:size(vChunk, 1)

        px = vChunk(v, 1);
        py = vChunk(v, 2);
        pz = vChunk(v, 3);
        lbl = lChunk(v);

        if lbl == 0, continue; end  % skip unlabelled (medial wall etc.)

        % find voxels within radius
        xi = find(abs(xv - px) <= radius);
        yi = find(abs(yv - py) <= radius);
        zi = find(abs(zv - pz) <= radius);

        if isempty(xi) || isempty(yi) || isempty(zi), continue; end

        % compute distances and Gaussian weights for nearby voxels
        [gx, gy, gz] = ndgrid(xv(xi), yv(yi), zv(zi));
        dist2 = (gx - px).^2 + (gy - py).^2 + (gz - pz).^2;
        w = exp(-dist2 / (2 * sigma^2));

        % accumulate weights for this parcel label
        voxWeights(xi, yi, zi, lbl) = voxWeights(xi, yi, zi, lbl) + w;

    end

    if mod(vStart, 50000) == 1
        fprintf('  ...processed %d / %d vertices\n', vEnd, nVerts);
    end

end

%% --- Assign each voxel the parcel with max accumulated weight ---
fprintf('Assigning parcel labels by majority vote...\n');
[maxW, volParc] = max(voxWeights, [], 4);   % nX x nY x nZ

% Zero out voxels with no votes at all
volParc(maxW == 0) = 0;

volParc = uint16(volParc);

fprintf('Done. Unique parcel labels assigned: %d\n', numel(unique(volParc)) - 1);
end