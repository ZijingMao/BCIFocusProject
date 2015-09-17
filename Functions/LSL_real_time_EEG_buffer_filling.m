function[EEG_buffer, op_status] =  ...
    LSL_real_time_EEG_buffer_filling(sampling_rate, EEG_buffer_size_sec, numchans)
% LSL_REAL_TIME_EEG_BUFFER_FILLING: This function is intended to be
% executed for a 'thread' by assign it as the callback function of a timer
% object. This function searches for a EEG stream of data (using LSL) adn
% constantly retrieves samples to fill a EEG buffer. Once the buffer is
% filled new samples replace old samples. Buffer filling time is also
% captured and save as global variables to be retrieved by the main script.

% Find EEG data stream

try
    disp('Looking for an EEG Stream...');
    lsl_object = lsl_loadlib();
    lsl_stream_handle = lsl_resolve_byprop(lsl_object, ...
        'name', 'Cog Cap-12 810');
    inlet_obj = lsl_inlet(lsl_stream_handle{1});
    [testing_sample, ~] = inlet_obj.pull_sample();
    disp('EEG Stream found. Testing sample collected');
    
catch
    
    % Notify error and exit
    op_status = false;
    EEG_buffer = [];
    disp('Stream Error: Verify BioSemi App and Data Streamming');
    disp('Exiting now...');
    return
end

%--- Here continues if EEG Data streamming was verified
% Pre-allocate data buffer
channels_to_extract = length(testing_sample);
EEG_buffer = zeros(channels_to_extract, (sampling_rate * EEG_buffer_size_sec));
buffer_sample_pos = 0;
buffer_number_samples = size(EEG_buffer, 2);

% Fill the EEG buffer
try
    tic
    for id=1:buffer_number_samples
        
        % Extract a sample
        [current_sample, ~] = inlet_obj.pull_sample();
        
        % Verify sample (empty or NaN raises flag to exit)
        if isempty(current_sample) || sum(isnan(current_sample)) > 0
            op_status = false;
            EEG_buffer = [];
            disp('Stream Error: Verify BioSemi App and Data Streamming');
            disp('Exiting now...');
            return;
        end
        
        % Verified sample now fills the EEG buffer
        buffer_sample_pos = buffer_sample_pos +1;
        EEG_buffer(:, buffer_sample_pos) = ...
            double(current_sample)';
        
    end
    
    % Remove extra channels and clean possible NaN's
    additional_chan = 0;  % Starts from 3, includes EX1 and EX2
%     EX1_pos = 34;         % This locations are based on the LSL output
%     EX2_pos = 35;         % when 32 channels are connected
    
    channel_positions = [1:(numchans + additional_chan)];
    % EEG_buffer = EEG_buffer(channel_positions, :);
    EEG_buffer( isnan(EEG_buffer) ) = 0;
    
    % Buffer should be filled, capture time and reset position counter
    EEG_buffer_fill_time = toc;
    disp(['Buffer filled in ', num2str(EEG_buffer_fill_time, 4), ' seconds']);
    op_status = true;
    return;
    
catch
    op_status = false;
    EEG_buffer = [];
    disp('Stream Error: Verify BioSemi App and Data Streamming');
    disp('Exiting now...');
end
end