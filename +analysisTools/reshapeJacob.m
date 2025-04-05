function jacob = reshapeJacob(jacob)

% Reshapes Jacobian matrix if necessary

    if length(size(jacob.A))>2  
        [nWL, nMeas, nVox]=size(jacob.A);
        jacob.A=reshape(permute(jacob.A,[2,1,3]),Nwl*nMeas,nVox);
    end

end