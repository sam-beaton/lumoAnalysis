% [s,tRange, sCh] = adaptedEnStimRejection(t,s,tIncAuto,tIncMan,tRange,tIncCh)
%
% Excludes stims that fall within the time points identified as 
% motion artifacts from HRF calculation.
%
%
% INPUT:
% t:        the time vector (#time points x 1)
% s:        s matrix (#time points x #conditions) containing 1 for 
%           each time point and condition that has a stimulus and 
%           zeros otherwise.
% tIncAuto: time points (#time points x 1) identified as motion 
%           artifacts by processing stream.
% tIncMan:  time points (#time points x 1) identified as motion 
%           artifacts by user.
% tRange:   an array of 2 numbers [t1 t2] specifying how many 
%           seconds surrounding motion artifacts, tIncMan and tIncAuto, 
%           to consider as excluded data and therefore exclude any stimuli 
%           which fall within those buffers.
%           Typically values are t1=-2 and t2 equal to the stimulus
%           duration.
% tIncCh:   a matrix with #time points x #channels, with 1's indicating data
%           included and 0's indicating motion artifacts on a channel by
%           channel basis              
%
% OUTPUT:
% s:        s matrix (#time points x #conditions) containing 1 for 
%           each time point and condition that has a stimulus that is 
%           included in the HRF calculation, -2 for a stimulus that is 
%           excluded automatically in the processing stream, -1 
%           for each stimulus excluded by a manually set patch and 
%           zeros otherwise.
% tRange:   same tRange array as in the input
% sCh:      s matrix (#time points x #channels) containing condition number
%           for each time point and condition that has a stimulus that is 
%           included in the HRF calculation in each channel, -2 for a 
%           stimulus that is excluded automatically in the processing 
%           stream, -1 for each stimulus excluded by a manually set patch 
%           and zeros otherwise.
%
% EDITED 30/3/25 by SLB:
% Takes channel-specific time inclusion and outputs channel specific stim
% inclusion

function [s, tRange, sCh] = adaptedEnStimRejection(t,s,tIncAuto,tIncMan,tRange,tIncCh)

dt = (t(end)-t(1))/length(t);
tRangeIdx = [floor(tRange(1)/dt):ceil(tRange(2)/dt)];

smax = max(s,[],2);
lstS = find(smax==1);

%initialise channel specific stim timings
[rowStim, colStim] = find(s == 1);
stimOccs = [rowStim, colStim];
stimOccs = sortrows(stimOccs, 1);

stimTimings = zeros(size(t));
for iStim = 1:length(lstS)
    stimTimings(lstS(iStim)) = stimOccs(iStim, 2);
end
sCh = repmat(stimTimings, 1, size(tIncCh, 2));



for iS = 1:length(lstS)
    lst = round(min(max(lstS(iS) + tRangeIdx,1),length(t)));
    if ~isempty(tIncAuto) && min(tIncAuto(lst))==0
        s(lstS(iS),:) = -2*abs(s(lstS(iS),:));
    end
    if ~isempty(tIncMan) && min(tIncMan(lst))==0
        s(lstS(iS),:) = -1*abs(s(lstS(iS),:));
    end
    if ~isempty(tIncCh)
        for iCh = 1:size(tIncCh, 2)
            if min(tIncCh(lst, iCh))==0
                sCh(lstS(iS), iCh) = -2;
            end
        end
    end
end
