function [BA_out, BSTD_out, BT_out, blocks] = getDataAverageBlock(data_in, pulse, dtPre, dtAfter, Tkeep)
% getDataAverageBlock Extracts stimulus-aligned data blocks including pre-stimulus baseline.
%
% INPUTS:
%   data_in - Input data array with time as the last dimension
%   pulse   - Vector of stimulus onset times (sample indices)
%   dtPre   - Number of samples to include before each pulse
%   dtAfter - Number of samples to include after each pulse
%   Tkeep   - Optional temporal mask; false values set corresponding timepoints to NaN
%
% OUTPUTS:
%   BA_out   - Block-averaged responses (demeaned across time)
%   BSTD_out - Standard deviation across blocks
%   BT_out   - Normalized block averages (BA_out ./ BSTD_out)
%   blocks   - Raw extracted blocks [channels x time x blocks]

%% Parameters and Initialization.
dims = size(data_in);
Nt = dims(end); % Assumes time is always the last dimension.
NDtf = (ndims(data_in) > 2); %#ok<ISMAT>
if ~exist('Tkeep','var'), Tkeep = true(Nt, 1); end
pulse = pulse(:);
Nbl = length(pulse);
dt = dtPre + dtAfter;

% Reshape for N-D input
if NDtf
    data_in = reshape(data_in, [], Nt);
end

% Apply temporal mask
data_in(:, ~Tkeep) = NaN;

% Allocate block matrix
Nm = size(data_in, 1);
blocks = NaN(Nm, dt, Nbl);

%% Extract blocks
for k = 1:Nbl
    startIdx = pulse(k) - dtPre;
    endIdx = pulse(k) + dtAfter - 1;

    % Adjust for pre-stimulus underflow
    padPre = 0;
    if startIdx < 1
        padPre = 1 - startIdx;
        startIdx = 1;
    end

    % Adjust for post-stimulus overflow
    padPost = 0;
    if endIdx > Nt
        padPost = endIdx - Nt;
        endIdx = Nt;
    end

    dataSegment = data_in(:, startIdx:endIdx);

    % Padding if needed
    if padPre > 0
        dataSegment = [NaN(Nm, padPre), dataSegment];
    end
    if padPost > 0
        dataSegment = [dataSegment, NaN(Nm, padPost)];
    end

    % Mark block invalid if any channel contains NaNs
    for iChan = 1:Nm
        if any(isnan(dataSegment(iChan, :)))
            dataSegment(iChan, :) = NaN;
        end
    end

    blocks(:, :, k) = dataSegment;
end

%% Compute block means and normalised output
BA_out = nanmean(blocks, 3);                    % Average across blocks
BSTD_out = nanstd(blocks, 0, 3);                % Std across blocks
BA_out = bsxfun(@minus, BA_out, nanmean(BA_out, 2));  % Demean each channel
BT_out = BA_out ./ BSTD_out;                   % Normalised
BT_out(~isfinite(BT_out)) = 0;

%% Reshape output to match input dimensions
if NDtf
    BA_out = reshape(BA_out, [dims(1:end-1), dt]);
    BSTD_out = reshape(BSTD_out, [dims(1:end-1), dt]);
    BT_out = reshape(BT_out, [dims(1:end-1), dt]);
    blocks = reshape(blocks, [dims(1:end-1), dt, Nbl]);
end
