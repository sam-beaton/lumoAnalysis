function parcelSensOut = getParcelSensitivity(maskParc, fooV, lightSensitivityMin, parcPercentMin)

% Returns wavelength-dependent parcel masks, containing those parcels which 
% the light coverage mask is sensitive to, for the provided coverage
% percentage.

    if ~exist('parcPercent', 'var')
        parcPercentMin = 0.5;
    elseif parcPercentMin >= 1 && parcPercentMin <= 100
        parcPercentMin = parcPercentMin./10;
    elseif parPercentMin >= 0 && parcPercentMin < 1
        %do nothing
    else 
        error('Please provide a valid spatial coverage percentage for your given parcellation')
    end

    % Find indices of parcellated GM
    indRegMaskParc = maskParc ~= 0;

    % Initialize output
    parcelSens = struct;
    
    for j = 1:numel(fieldnames(fooV))

        %find indices of light coverage mask
        lambda = fooV.(['lambda' num2str(j)]);
        sensMask = lambda > lightSensitivityMin;
        
        % Sensitivity mask > threshold, restricted to GM
        parcLight = sensMask & indRegMaskParc;
        
        % obtain GM space covered by array
        indRegParcBin  = zeros([size(maskParc)]);
        indRegParcBin(parcLight) = 1; %create binary mat. to multiply by FooV
        
        % parcels overlapping with sensitivity mask
        parcFooVOnly = indRegParcBin.*maskParc;
        parcelsIncl = unique(parcFooVOnly); %find all parcel numbers in a vector
        parcelsIncl(parcelsIncl == 0) = []; %remove 0 as it corresponds to non-GM space

        %initialise output mask
        parcelsSens = zeros(size(maskParc)); 

        for iParc = length(parcelsIncl):-1:1
            %check proportion of parcel voxels covered by array
            if length(find(parcFooVOnly == parcelsIncl(iParc))) / length(find(maskParc == parcelsIncl(iParc))) >= parcPercentMin
                %if >= 50%, include parcel in mask
                parcelsSens(maskParc == parcelsIncl(iParc)) = parcelsIncl(iParc);
            else
                %otherwise remove it from list of included parcels
                parcelsIncl(iParc)=[];
            end
        end
        if j == 1
            parcelSensOut.lambda1 = parcelsSens;
        elseif j == 2
            parcelSensOut.lambda2 = parcelsSens;
        end
    end
end
