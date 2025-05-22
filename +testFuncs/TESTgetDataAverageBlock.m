function [BA_out, BSTD_out, BT_out, blocks] = TESTgetDataAverageBlock(data_in, pulse, dt, Tkeep)

    % Dimensions
    dims = size(data_in);
    Nt = dims(end);
    isND = (ndims(data_in) > 2);
    if nargin < 4 || isempty(Tkeep), Tkeep = true(Nt, 1); end

    % Trim pulses that exceed bounds
    Nbl = sum((pulse + dt - 1) <= Nt);
    pulse = pulse(1:Nbl);

    % Reshape if needed
    if isND
        data_in = reshape(data_in, [], Nt);  % [channels x time]
    end

    % Apply temporal mask
    data_in(:, ~Tkeep) = NaN;

    % Pre-allocate
    [nCh, ~] = size(data_in);
    blocks = NaN(nCh, dt, Nbl);

    % Vectorized block extraction (loop over blocks only)
    for k = 1:Nbl
        idxEnd = pulse(k) + dt - 1;
        if idxEnd <= Nt
            segment = data_in(:, pulse(k):idxEnd);
        else
            segment = [data_in(:, pulse(k):end), NaN(nCh, idxEnd - Nt)];
        end

        % Set entire row to NaN if any NaNs exist in that row
        nanRows = any(isnan(segment), 2);
        segment(nanRows, :) = NaN;

        blocks(:, :, k) = segment;
    end

    % Compute block averages and standard deviations
    BA_out = nanmean(blocks, 3);
    BSTD_out = nanstd(blocks, 0, 3);

    % Demean and normalize
    BA_out = BA_out - nanmean(BA_out, 2);
    BT_out = BA_out ./ BSTD_out;
    BT_out(~isfinite(BT_out)) = 0;

    % Restore original dimensions if needed
    if isND
        outDims = [dims(1:end-1), dt];
        BA_out = reshape(BA_out, outDims);
        BSTD_out = reshape(BSTD_out, outDims);
        BT_out = reshape(BT_out, outDims);
        blocks = reshape(blocks, [dims(1:end-1), dt, Nbl]);
    end
end
