% This script computes the error rate based on CCA's frequency detection of
% a set of SSVEP stimulus. Once a SSVEP experiment has been recorded,
% configure this script to compute the experiment's performance. It can
% also test multiple datasets and multiples parameters using a parfor loop.

clear
close all
clc

% Add EEGlab and needed file\folders to the path
recordings_folder = 'Frequency_set_2_new_layout';  % Specific recording folder to use
use_parfor = false;                                % Activate or not parallel for loops
[dataset_folder, script_folder] = set_path_and_parpool(use_parfor, recordings_folder);

%% Define the default SSVEP Channel/label layout (11 Channels + Refs)
ssvep_EEG_base_10_20_position = [21:23, 25:32];
ssvep_EEG_LSL_matrix_position = 3:13;
ssvep_EEG_channel_10_20_position = [21:23, 25:32];
ssvep_EEG_channel_10_20_labels = {'Pz', 'POz', 'PO3', 'PO4', 'Oz', ...
    'O1', 'O2', 'P3', 'P4', 'PO7', 'PO8'};

%% Configure FIXED input parameters

% Configure file locations
input_configuration.dataset_folder = dataset_folder;
input_configuration.script_folder = script_folder;

% Configure data recording properties
input_configuration.common_reference = true; % Enable/disable (true/false) common referencing
input_configuration.reference_channel = [34, 35];   % Either 2(Cz) or [34, 35](EX1 and EX2)
input_configuration.sampling_rate = 512;     % Sampling rate when recording

% Configure data pre-procesing and sliding window
input_configuration.band_pass_filter_freqs = [0.3, 60]; % Band-pass filtering
input_configuration.step_size_time = 0.5;             % Sliding window step
input_configuration.event_data_extraction_offset_time = 1;
input_configuration.event_start_number = 2;

% Configure SSVEP stimulus properties
input_configuration.target_event_labels = {'Up', 'Down', 'Left', 'Right', 'Face'};
input_configuration.ssvep_freqs = [6, 7.5, 8.57, 10, 10.9]; % Stimulus freqs
input_configuration.target_event_freqs = [7.5, 8.57, 10, 6, 10.9];

%% Configure Sets of parameter values (each combination will be evaluated)
input_configuration.participant_filename = {'mauricio_new_layout.xdf', ...
    'tinghe_new_layout_2.xdf', 'zijing_new_layout.xdf'};
input_configuration.channel_subset = 11;
input_configuration.window_size_time = 1:0.5:4;
input_configuration.num_harmonics = 2:6;

%% Pre-allocate output error rates, compute performance

% Get the number of values for each parameter
n_participant = length(input_configuration.participant_filename);
n_window = length(input_configuration.window_size_time);
n_harmonics = length(input_configuration.num_harmonics);
n_channel_subset = length(input_configuration.channel_subset);

% Create error rate table template
error_rate_table = zeros((n_window*n_harmonics*n_channel_subset)+1, ...
    length(input_configuration.target_event_labels)+3);

% Fill the first three columns of the table, which are fixed values
error_rate_table(1, 1) = n_harmonics;
error_rate_table(1, 2) = n_window;
error_rate_table(1, 3) = n_channel_subset;
error_rate_table(1, 4:end) = input_configuration.target_event_freqs;

% Create also tables to store the computed number of data events/epochs
data_events_table = error_rate_table;
frequency_epochs_Table = error_rate_table;

% This cell will contain a table for each participant
participant_freq_error_rate = cell(1, length(input_configuration.participant_filename));
shadow_conf_struct = input_configuration; % Struct copy to input to fcn.

% Loop through configurations, for each case compute and save error rate
for loop_participant=1:n_participant
    
    % Reset the performance table for current participant
    participant_freq_error_rate{loop_participant} = error_rate_table;
    row_counter = 1;
    
    % Replace inputs on shawdow config. struct
    shadow_conf_struct.participant_filename = ...
        input_configuration.participant_filename(loop_participant);
    
    for loop_channel=1:n_channel_subset
        
        % Replace inputs on shadow config. struct
        % Replace number of channels subset
        shadow_conf_struct.channel_subset = ...
            input_configuration.channel_subset(loop_channel);
        
        % Replace LSL matrix channel positions (current subset)
        shadow_conf_struct.EEG_LSL_matrix_position = ...
            ssvep_EEG_LSL_matrix_position(1:shadow_conf_struct.channel_subset);
        
        % Replace 10-20 channel locations (current subset)
        shadow_conf_struct.EEG_channel_10_20_position = ...
            ssvep_EEG_channel_10_20_position(1:shadow_conf_struct.channel_subset);
        
        % Replace channel labels (current subset)
        shadow_conf_struct.EEG_channel_10_20_labels = ...
            ssvep_EEG_channel_10_20_labels(1:shadow_conf_struct.channel_subset);
        
        for loop_harmonic=1:n_harmonics
            
            % Replace inputs on shawdow config. struct
            shadow_conf_struct.num_harmonics = ...
                input_configuration.num_harmonics(loop_harmonic);
            
            for loop_window=1:n_window
                
                % Replace inputs on shawdow config. struct
                shadow_conf_struct.window_size_time = ...
                    input_configuration.window_size_time(loop_window);
                
                %----------------------------------------------------------
                % Compute performance, save values on next row on table
                [stimulus_error_rates, event_counter, event_num_epochs] = ...
                    compute_ssvep_dataset_cca_performance(shadow_conf_struct);
                %----------------------------------------------------------
                
                % Update table
                % Update row counter
                row_counter = row_counter + 1;
                
                % Update current number of CCA harmonics
                participant_freq_error_rate{loop_participant}(row_counter, 1) = ...
                    input_configuration.num_harmonics(loop_harmonic);
                
                % Update current sliding window length
                participant_freq_error_rate{loop_participant}(row_counter, 2) = ...
                    input_configuration.window_size_time(loop_window);
                
                % Update current number of channels subset
                participant_freq_error_rate{loop_participant}(row_counter, 3) = ...
                    input_configuration.channel_subset(loop_channel);
                
                % Update error rates
                participant_freq_error_rate{loop_participant}(row_counter, 4:end) = ...
                    stimulus_error_rates;
                
                % -------------- Update events/epoch tables ---------------
                data_events_table(row_counter, 4:end) = event_counter;
                frequency_epochs_Table(row_counter, 4:end) = event_num_epochs;
                
                
            end
        end
    end
end

disp('Process completed!');


