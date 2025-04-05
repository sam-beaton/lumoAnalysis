function fooV = getLightCoverage(A, iA, infoA, infoData, params)

 
% Takes a Jacobian matrix and returns the light coverage 

    nLambda = length(unique(infoData.pairs.lambda)); %number of wavelengths
    if ~exist('dim', 'var')
        dim = infoA.tissue.dim;
    end

    % invert A
    for j = 1:nLambda
        keep = (infoA.pairs.WL == j) & (infoA.pairs.r3d <= params.maxChannelDistance) & infoData.MEAS.GI;
        aWav = squeeze(A(keep,:));
        if j == 1
            ffr=makeFlatFieldRecon(aWav, iA.lambda1);
            fooV1=Good_Vox2vol(ffr,dim);
            fooV.lambda1=fooV1./max(fooV1(:));
        elseif j == 2
            ffr=makeFlatFieldRecon(aWav, iA.lambda2);
            fooV2=Good_Vox2vol(ffr,dim);
            fooV.lambda2=fooV2./max(fooV2(:));
        else
            error('Invalid value for j: %d. Expected 1 or 2.', j);
        end
    end
end