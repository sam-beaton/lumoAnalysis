function fooV = getLightCoverage(A, iJacob, info, params)

% Takes a Jacobian matrix and returns the light coverage 

    nLambda = length(unique(info.pairs.lambda)); %number of wavelengths
    if ~exist('dim', 'var')
        dim = info.tissue.dim;
    end

    % invert A
    for j = 1:nLambda
        keep = (info.pairs.WL == j) & (info.pairs.r3d <= params.maxChannelDistance) & info.MEAS.GI;
        aWav = squeeze(A(keep,:));
        if j = 1
            ffr=makeFlatFieldRecon(aWav,iJacob.lambda1);
            fooV1=Good_Vox2vol(ffr,dim);
            fooV.lambda1=fooV1./max(fooV1(:));
        elseif j = 2
            ffr=makeFlatFieldRecon(aWav,iJacob.lambda2);
            fooV2=Good_Vox2vol(ffr,dim);
            fooV.lambda2=fooV2./max(fooV2(:));
        else
            error('Invalid value for j: %d. Expected 1 or 2.', j);
        end
    end
end