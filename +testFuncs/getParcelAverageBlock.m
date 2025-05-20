function [BA_out, BSTD_out, BT_out, blocks] = getParcelAverageBlock(data_in, pulse, dt, Tkeep)

% BLOCKAVERAGE Averages data by stimulus blocks.
%
%   data_out = BLOCKAVERAGE(data_in, pulse, dt) takes a data array
%   "data_in" and uses the pulse and dt information to cut that data 
%   timewise into blocks of equal length (dt), which are then averaged 
%   together and output as "data_out".
% 
%   Tkeep is a temporal mask. Any time points with a zero in this vector is
%   set to NaN.
%
% Copyright (c) 2017 Washington University 
% Created By: Adam T. Eggebrecht
% Eggebrecht et al., 2014, Nature Photonics; Zeff et al., 2007, PNAS.
%
% Washington University hereby grants to you a non-transferable, 
% non-exclusive, royalty-free, non-commercial, research license to use 
% and copy the computer code that is provided here (the Software).  
% You agree to include this license and the above copyright notice in 
% all copies of the Software.  The Software may not be distributed, 
% shared, or transferred to any third party.  This license does not 
% grant any rights or licenses to any other patents, copyrights, or 
% other forms of intellectual property owned or controlled by Washington 
% University.
% 
% YOU AGREE THAT THE SOFTWARE PROVIDED HEREUNDER IS EXPERIMENTAL AND IS 
% PROVIDED AS IS, WITHOUT ANY WARRANTY OF ANY KIND, EXPRESSED OR 
% IMPLIED, INCLUDING WITHOUT LIMITATION WARRANTIES OF MERCHANTABILITY 
% OR FITNESS FOR ANY PARTICULAR PURPOSE, OR NON-INFRINGEMENT OF ANY 
% THIRD-PARTY PATENT, COPYRIGHT, OR ANY OTHER THIRD-PARTY RIGHT.  
% IN NO EVENT SHALL THE CREATORS OF THE SOFTWARE OR WASHINGTON 
% UNIVERSITY BE LIABLE FOR ANY DIRECT, INDIRECT, SPECIAL, OR 
% CONSEQUENTIAL DAMAGES ARISING OUT OF OR IN ANY WAY CONNECTED WITH 
% THE SOFTWARE, THE USE OF THE SOFTWARE, OR THIS AGREEMENT, WHETHER 
% IN BREACH OF CONTRACT, TORT OR OTHERWISE, EVEN IF SUCH PARTY IS 
% ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
%
% See Also: NORMALIZE2RANGE_TTS.

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
blocks=zeros(Nm,dt,Nbl);
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