function Tkeep = findtKeep(data_in, params, info)

% Finds the timepoints which have not ben excluded due to motion and which
% should therefore be included in analysis.

%% Parameters and Initialization.
dims = size(data_in); %'data': (channels)*(sample points)
Nt = dims(end); % last sample point - Assumes time is always the last dimension.
NDtf = (ndims(data_in) > 2); % sets to one if data >2 dimensional (needs transforming)

% account for different sampling frequencies
dtPre = floor(params.dtPre*(info.system.framerate/10)); 
dtAfter = floor(params.dtAfter*(info.system.framerate/10));

%% N-D Input (for 3-D or N-D voxel spaces).
if NDtf
    data_in = reshape(data_in, [], Nt);
end

%% only use default standard temporal mask if not provided and
% info.paradigmFull.sCh doesn't exist
if ~exist('Tkeep','var')
    Tkeep=ones(size(data_in)); 
    if isfield(info, 'paradigmFull') && isfield(info.paradigmFull, 'sCh')

        % find non-baseline stim markers
        pulse = [1:length(info.paradigmFull.synchpts)]';
        pulse = setdiff(pulse, info.paradigmFull.Pulse_1);


        %find stim markers where motion was detected (negative values)
        [samplepoints, channels] = find(info.paradigmFull.sCh < 0);
        negStims = [samplepoints, channels];
        negStims = sortrows(negStims);

        %select negative stim markers corresponding to non-baseline 
        % pulse/stimulus
        pulseRows = ismember(negStims(:,1), info.paradigmFull.synchpts(pulse));
        negLocs = negStims(pulseRows, :);
        
        % set blocks in Tkeep to 0 if a negative stim number
        for i = 1:size(negLocs, 1)
            stimTime = negLocs(i, 1);  % stimulus time
            stimChan = negLocs(i, 2);  % channel
            % Set start & ensure index is within bounds
            startRow = max(stimTime, 1);   
            % Set end & ensure index is within bounds
            endRow = min(stimTime + dtAfter, size(data_in, 2)); 
            % Set the specified column (channel) to zero in the selected range
            Tkeep(stimChan, startRow:endRow) = 0;
        end
    end
end

end


