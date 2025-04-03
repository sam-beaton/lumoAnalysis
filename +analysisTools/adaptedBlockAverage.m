function [BA_out,BSTD_out,BT_out,blocks] = adaptedBlockAverage(data_in, params, info, Tkeep)

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
%
% Edited 1/4/25 SLB:
% Accounts for trial rejection on a channel-by-channel basis using
% info.paradigmFull.sCh (see 'getNdotFile.m')
% Incorporates a baseline into the block (for later use when testing
% statistical significance) and accounts for possibility this baseline may
% protrude beyond start of data in file.

%% Parameters and Initialization.
dims = size(data_in); %'data': (channels)*(sample points)
Nt = dims(end); % last sample point - Assumes time is always the last dimension.
NDtf = (ndims(data_in) > 2); % sets to one if data >2 dimensional (needs transforming)
pulse = [info.paradigmFull.Pulse_2; info.paradigmFull.Pulse_3; info.paradigmFull.Pulse_4];
Nbl = length(pulse); % number of presentations for first stim - based on full paradigm

% Check to make sure that the block after the last synch point for this
% pulse does not exceed the data's time dimension. 
if (info.paradigmFull.synchpts(pulse(end)) + params.dtAfter) > Nt
    pulse = pulse(1:end-1);
    Nbl = Nbl - 1; 
end
% Check to make sure that the block before the first synch point for this
% pulse does not exceed the data's time dimension. 
if (info.paradigmFull.synchpts(pulse(1)) - params.dtPre) < 1
    pulse = pulse(2:end);
    Nbl = Nbl - 1; 
end

%% N-D Input (for 3-D or N-D voxel spaces).
if NDtf
    data_in = reshape(data_in, [], Nt);
end

%% only use default standard temporal mask if not provided and
% info.paradigmFull.sCh doesn't exist
if ~exist('Tkeep','var')
    Tkeep=ones(size(data_in)); 
    if isfield(info, 'paradigmFull') && isfield(info.paradigmFull, 'sCh')
        %find stim markers where motion was detected (negative values)
        [rows, cols] = find(info.paradigmFull.sCh < 0);
        negStims = [rows, cols];
        negStims = sortrows(negStims);

        %select negative stim markers corresponding to pulse/stimulus for
        %averaging
        pulseRows = ismember(negStims(:,1), info.paradigmFull.synchpts(pulse));
        pulseLocs = negStims(pulseRows, :);
        
        % set blocks in Tkeep to 0 if a negative stim number
        for i = 1:size(pulseLocs, 1)
            stimTime = pulseLocs(i, 1);  % stimulus time
            stimChan = pulseLocs(i, 2);  % channel
            startRow = max(stimTime - params.dtPre, 1);   % Ensure index is within bounds
            endRow = min(stimTime + params.dtAfter, size(data_in, 1)); % Ensure index is within bounds
            % Set the specified column (channel) to zero in the selected range
            Tkeep(stimChan, startRow:endRow) = 0;
        end
    end
end

%% Incorporate Tkeep  
data_in(~Tkeep)=NaN; %sets unwanted samplepoints to NaN 


%% Cut data into blocks.
Nm=size(data_in,1); %number of channels
dt = params.dtPre+params.dtAfter; % number of samples in block
blocks=zeros(Nm,dt,Nbl); % (numChans)*(sample points post-stim)*(number of stim presentations)
for k = 1:Nbl
    synchSamp = info.paradigmFull.synchpts(pulse(k)); % timing of pulse/stim
    blockStart = synchSamp - params.dtPre; % Start index of the block
    blockEnd = synchSamp + params.dtAfter - 1; % Block end index
    dtbPre = 0; 
    dtbPost = 0;

    % if block extends before start of the recording:
    if blockStart < 1
        dtbPre = abs(blockStart) + 1; % get # missing time points @ start
        blockStart = 1; % adjust 'start' to file start
    else
        dtbPre = 0; 
    end

    % if block extends beyond end of the recording
    if blockEnd > Nt
        dtbPost = blockEnd - Nt; % get # missing time points @ end
        blockEnd = Nt; % adjust 'end' to file end
    else
        dtbPost = 0; 
    end

    % extract available data
    dataSegment = data_in(:, blockStart:blockEnd); 

    % add NaN padding where needed via concatenation 
    blocks(:, :, k) = cat(2, NaN(Nm, dtbPre), dataSegment, NaN(Nm, dtbPost)); 

    % Original structure, not finished (only accounts for block extending
    % beyond end of file):

%     if (synchSamp + dt - 1) <= Nt %checks block contained within recording
%         blocks(:, :, k) = data_in(:, synchSamp-params.dtPre:synchSamp+params.dtAfter - 1); 
%     else
%         dtb = synchSamp + dt - 1 - Nt; % ADAPT
%         blocks(:, :, k) = cat(2,data_in(:,synchSamp:end),NaN(size(data_in,1),dtb)); % ADAPT
%     end
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