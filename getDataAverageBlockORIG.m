function [BA_out, BSTD_out, BT_out, blocks] = getDataAverageBlock(data_in, pulse, dt, Tkeep)

% getDataAverageBlock Averages data by stimulus blocks.
%
% GETPARCELAVERAGEBLOCK Extracts stimulus-aligned data blocks and computes statistics.
%
%   Segments time-series data into stimulus-aligned blocks of fixed duration,
%   computes block-averaged responses, and returns demeaned/normalised outputs.
%
% INPUTS:
%   data_in - Input data array with time as the last dimension
%   pulse   - Vector of stimulus onset times (sample indices)
%   dt      - Block duration in samples
%   Tkeep   - Optional temporal mask; false values set corresponding timepoints to NaN
%
% OUTPUTS:
%   BA_out   - Block-averaged responses (demeaned across time)
%   BSTD_out - Standard deviation across blocks
%   BT_out   - Normalized block averages (BA_out ./ BSTD_out)
%   blocks   - Raw extracted blocks [channels x time x blocks]
%
% NB:
%   - Time dimension is assumed to be the last dimension of data_in
%   - Blocks containing any NaN values are excluded from averaging
%   - Supports both 2D (channels x time) and N-D input arrays
%   - Pulses extending beyond data boundaries are handled with NaN padding
%
% SLB 22/5/25
% Edited from BLOCKAVERAGE in the NeuroDOT toolbox

%% Parameters and Initialization.
dims = size(data_in);
Nt = dims(end); % Assumes time is always the last dimension.
NDtf = (ndims(data_in) > 2); %#ok<ISMAT>
Nbl = length(pulse);
if ~exist('Tkeep','var'), Tkeep=ones(Nt,1)==1;end % temporal mask

% Check to make sure that the block after the last synch point for this
% pulse does not exceed the data's time dimension. 
if (dt + pulse(end)-1) > Nt
    Nbl = Nbl - 1;
end

%% N-D Input (for 3-D or N-D voxel spaces).
if NDtf
    data_in = reshape(data_in, [], Nt);
end

%% Incorporate Tkeep  
data_in(:,~Tkeep)=NaN;


%% Cut data into blocks.
Nm=size(data_in,1);
blocks = NaN(Nm, dt, Nbl);
for k = 1:Nbl
    if (pulse(k) + dt - 1) <= Nt
        % extract the data segment
        dataSegment = data_in(:, pulse(k):(pulse(k) + dt - 1));
        
        % check if the segment contains any NaNs (channel-wise)
        for iChan = 1:size(dataSegment, 1)
            if any(isnan(dataSegment(iChan, :)))
                % if so, set the entire channel to NaN
                dataSegment(iChan, :) = NaN;
            else
                %dataSegment(iChan, :) = dataSegment(iChan, :) - nanmean(dataSegment(iChan, :)); 
            end
        end
        % assign
        blocks(:, :, k) = dataSegment;
    else % For blocks extending beyond data end:
        
        % number of sampled values to be padded with NaNs
        dtb = pulse(k) + dt - 1 - Nt;
        
        % extrac available data
        dataSegment = data_in(:, pulse(k):end);
        
        % check if the segment contains any NaNs (channel-wise)
        for iChan = 1:size(dataSegment, 1)
            if any(isnan(dataSegment(iChan, :)))
                % if so, set the entire channel to NaN
                dataSegment(iChan, :) = NaN;
            else
                %dataSegment(iChan, :) = dataSegment(iChan, :) - nanmean(dataSegment(iChan, :)); 
            end
        end
        
        % concatenate with NaN padding
        blocks(:, :, k) = cat(2, dataSegment, NaN(size(data_in, 1), dtb));
    end
end

%% Average blocks and return.
BA_out = nanmean(blocks, 3);
BSTD_out = nanstd(blocks, [],3);
BA_out = bsxfun(@minus,BA_out,nanmean(BA_out,2));
BT_out = BA_out./BSTD_out;
BT_out(~isfinite(BT_out))=0;


%% N-D Output.
if NDtf
    BA_out = reshape(BA_out, [dims(1:end-1), dt]);
    BSTD_out = reshape(BSTD_out, [dims(1:end-1), dt]);
    BT_out = reshape(BT_out, [dims(1:end-1), dt]);
    blocks=reshape(blocks,[dims(1:end-1), dt,Nbl]);
end