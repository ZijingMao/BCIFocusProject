function [event_error_rate, event_counter, event_num_epochs] = ...
    compute_ssvep_dataset_cca_performance(input_configuration)
%SSVEP_DATASET_CCA_PERFORMANCE Extracts data chunks for input events on
%input datasets and then epoch each chunk using a sliding window. For each
%epoc CCA detection is applied and incorrect  decisions are marked to
%compute error rate among all epochs.

% Extract variables from the configuration struct
participant_filename = input_configuration.participant_filename;
sampling_rate = input_configuration.sampling_rate;
band_pass_filter_freqs = input_configuration.band_pass_filter_freqs;
step_size_time = input_configuration.step_size_time;
window_size_time = input_configuration.window_size_time;
num_harmonics = input_configuration.num_harmonics;
target_event_labels = input_configuration.target_event_labels;
ssvep_freqs = input_configuration.ssvep_freqs;
target_event_freqs = input_configuration.target_event_freqs;
event_data_extraction_offset_time = input_configuration.event_data_extraction_offset_time;
event_start_number = input_configuration.event_start_number;
dataset_folder = input_configuration.dataset_folder;
reference_channel = input_configuration.reference_channel;
common_reference = input_configuration.common_reference;
EEG_channel_10_20_position = input_configuration.EEG_channel_10_20_position;
EEG_LSL_matrix_position = input_configuration.EEG_LSL_matrix_position;


% Check for extraction of ALL participants (current data sub-folder)
n_input_participants = length(participant_filename);
filename_queque = participant_filename;

if n_input_participants == 1 && strcmp(participant_filename{1}, 'all')
    
    % Here all participant were requested, create a list of files
    cd(dataset_folder);
    data_folder_filelist = dir(pwd);
    data_folder_filelist = extractfield(data_folder_filelist, 'name');
    dataset_counter = 0;
    
    % Find valid .XDF files
    for loop=1:length(data_folder_filelist)
        
        format_check = strfind(data_folder_filelist{loop}, '.xdf');
        format_check = isnumeric(format_check);
        
        if format_check > 1
            % Increase counter and extract the filename
            dataset_counter = dataset_counter +1;
            filename_queque{dataset_counter} = data_folder_filelist{loop};
        end
    end
    
    % Update the number of participants
    n_input_participants = dataset_counter;
end

% Pre-allocate and extract data from .XDF files
participant_data = cell(1, n_input_participants);
participant_num_events = zeros(1, n_input_participants);

for loop=1:n_input_participants
    
    % Load current filename and extract data
    participant_data{loop} = pop_loadxdf(participant_filename{loop}, ...
        'streamtype', 'EEG', 'exclude_markerstreams', {});
    
    % Extract current participant number of events
    participant_num_events(loop) = length(participant_data{loop}.event);
end

% Define and pre-allocate variables
num_event_targets = length(target_event_labels);
n_events_all = sum(participant_num_events);
event_counter = zeros(1, num_event_targets);
flicker_events_data = cell(num_event_targets, n_events_all);

% Extract events data by target
for id1=1:n_input_participants
    
    % Compute extraction offset samples
    event_data_extraction_offset_points = ...
        participant_data{id1}.srate * event_data_extraction_offset_time;
    
    
    % Extract all events data for current participant
    current_data = participant_data{id1}.data;
    current_events = participant_data{id1}.event;
    
    % Loop through events
    for id2=event_start_number:length(participant_data{1, id1}.event) - 1
        
        % Extract current and next event labels and time point position
        current_event_label = current_events(id2).type;
        next_event_label = current_events(id2+1).type;
        current_event_npoint = floor(current_events(id2).latency);
        next_event_npoint = floor(current_events(id2+1).latency);
        
        % Check for any of the target events
        for id3=1:num_event_targets
            
            if strcmp(current_event_label, target_event_labels{id3}) && ...
                    strcmp(next_event_label, 'Cross')
                
                % Target Flickering event found, increase counter and
                % extract data
                event_counter(id3) = event_counter(id3) + 1;
                flicker_events_data{id3, event_counter(id3)} = double(...
                    current_data(:, ...
                    (current_event_npoint+event_data_extraction_offset_points):...
                    next_event_npoint));
                
            end
        end
    end
end

%% Perform sliding window on each event data chunk

% Define variables
step_size_points = sampling_rate * step_size_time;
window_size_points = sampling_rate * window_size_time;
event_epoched_data = cell(1, num_event_targets);
event_num_epochs = zeros(1, num_event_targets);
num_slides = zeros(1, num_event_targets);

% Loop target events
for id1=1:num_event_targets
    
    % Compute the minimum size (time points) among all events
    [num_chan_array, num_points_array] = cellfun(@size, flicker_events_data(id1, 1:event_counter(id1)));
    min_size = min(num_points_array);
    
    % Get number of slides and pre-allcated windowed event data (epochs)
    num_slides(id1) = floor((min_size - window_size_points) / step_size_points) + 1;
    event_num_epochs(id1) = event_counter(id1) * num_slides(id1);
    
    % Pre-allocate windowed epoch data (for current target event)
    event_epoched_data{id1} = ...
        zeros(min(num_chan_array), window_size_points, event_counter(id1), num_slides(id1));
    
    %--- Extract epochs (using sliding window)
    for id2=1:event_counter(id1)
        
        % Extract current event data chunk
        current_data = flicker_events_data{id1, id2};
        
        % Loop through window slides (as many as computed before)
        window_start_point = 0;
        for id3=1:num_slides(id1)
            
            % Extract windowed slides (epochs)
            window_start_point = window_start_point + 1;
            current_epoch = current_data(:, window_start_point:(window_start_point + window_size_points)-1);
            window_start_point = (window_start_point + step_size_points) - 1;
            
            % Save current epoch into the event epochs matrix
            event_epoched_data{id1}(:, :, id2, id3) = current_epoch;
            
        end
    end
    
end

%% Perform pre-processing and CCA detection on every epoch

% Pre-allocate variables
event_error_rate = zeros(1, num_event_targets);

% Loop through target events
for id1=1:num_event_targets
    
    % Pre-allocate predicted labels for current target event
    CCA_event_predicted_labels = zeros(1, event_num_epochs(id1));
    event_epoch_counter = 0;
    
    % Loop through every epoch (every "slide/window" per event data chunk)
    for id2=1:event_counter(id1)
        for id3=1:num_slides(id1)
            
            % Extract current epoch (it is redundant but the code is clear)
            current_epoch = event_epoched_data{id1}(:, :, id2, id3);
            event_epoch_counter = event_epoch_counter + 1;
            
            % Append reference channels (EX1, EX2 as the last two chans)
            current_epoch(end-1:end, :) = current_epoch(reference_channel, :);
            
            % Pre-process (filter, refer, etc.) EEG buffer
            [~, processed_epoch] = ...
                ssvepLiveProcessingAndDetection_v3(current_epoch, ...
                sampling_rate, band_pass_filter_freqs, ...
                common_reference, EEG_channel_10_20_position, ...
                EEG_LSL_matrix_position);
            
            % Perform CCA using the SSVEP stimulus frequencies
            [~, CCA_predicted_freq] = ...
                ssvep_cca_correlation(processed_epoch, ssvep_freqs, num_harmonics);
            
            if CCA_predicted_freq ~= target_event_freqs(id1)
                
                % Mark epoch as error epoch
                CCA_event_predicted_labels(event_epoch_counter) = 1;
            end
        end
    end
    
    % After peform prediction on all epochs, compute error rate
    event_error_rate(id1) = sum(CCA_event_predicted_labels)/ event_num_epochs(id1);
end

% Clear command window
clc
end
