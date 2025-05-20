function [cortexMuA, iJacobOut]  = imageReconstruction(data, jacobian, info, params)

% lmdata, jacob, data.info, paramsFile
% Performs image reconstruction, using 'reconstruct_img' from NeuroDOT
% toolbox

    nVox=size(jacobian.A,2); %number of voxels in voxellated sensitivity map
    nT = size(data,2); %number of time points in (resampled) recording
    nLambda = length(unique(info.pairs.lambda)); %number of wavelengths
    % initialise absorption matrix i.e. inversion of A, twice (once for each WL)
    cortexMuA = zeros(nVox, nT, nLambda);

    % invert A
    for j = 1:nLambda
        keep = (jacobian.info.pairs.WL == j) & (jacobian.info.pairs.r3d <= params.maxChannelDistance) & info.MEAS.GI; 
        %fprintf('Inverting Jacobian: wavelength %g\n', j)             
        iJacob = Tikhonov_invert_Amat(jacobian.A(keep, :), 0.1, 0.1); % Invert A-Matrix
        %fprintf('Smoothing Inverse: wavelength %g\n', j)   
        iJacob = analysisTools.adaptedSmoothAmat(iJacob, jacobian.info.tissue.dim, 3); % Smooth Inverted A-Matrix
        %fprintf('Reconstructing Volume: wavelength %g\n', j) 
        cortexMuA(:, :, j) = reconstruct_img(data(keep, :), iJacob); % Reconstruct Image Volume
        if j == 1
            iJacobOut.lambda1 = iJacob;
        else
            iJacobOut.lambda2 = iJacob;
        end
    end

end