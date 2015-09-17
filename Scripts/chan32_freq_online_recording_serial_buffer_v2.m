% Version 4, updated help will be added later
clear
close all
clc

% Add EEGlab and needed file\folders to the path
exe_folder = 'E:\Zijing\Documents\EEG\Cognionics\Realtime\';  % Specific recording folder to use
cd([exe_folder 'Scripts']);
use_parfor = false;                      % active or not parallel for loops
[dataset_folder, script_folder] = set_path_and_parpool(use_parfor, exe_folder);

%% SSVEP frequency detection fine-tunning
CCA_num_harmonics = 2;
EEG_buffer_size_sec = 3;          % EEG buffer Length (adjust per participant)
freq_range_spectra_vis = [1, 40]; % Range for Freq. spectrum visualization
band_pass_filter_freqs = [0.3, 50]; % EEGLAB- bandpass filtering

%% Open an Outlet for LSL streamming of marker events

% Use library an open outlet
lib = lsl_loadlib();
info_ = lsl_streaminfo(lib,'MyMarkerStream','Markers',1,0,'cf_string','MATLABProcessingEvents');
outlet = lsl_outlet(info_);

% Send event label (using LSL)
event_str = 'MATLAB_start';
outlet.push_sample({event_str});
pause(1);

%% Configure Channel locations, re-referincg and SSVEP frequencies
new_ssvep_channel_config = 0; % switch here between NEW, OLD or Subset SSVEP Channel layout
common_reference = true;
channel_length = 32;
EEG_channel_10_20_position = 1:channel_length;
EEG_LSL_matrix_position = 1:channel_length;
EEG_channel_10_20_labels = cell(1, channel_length);

% Flickering freq. on SSVEP presentation
% freq_band = 8:17;
% resolution = 1;
% freq_bank = freq_band(1):resolution:freq_band(2);
% xtick = freq_bank;

%% Launch real-time raw data plot (vis_stream)

%--- Configure and launch EEG buffer filling function, get filled buffer
% load EEG_channel_locations_32_10_20_system
sampling_rate = 500;              % Fron BioSemi device OR simulation
number_of_chanels = length(EEG_channel_10_20_labels); % Update this to match experiment

%--- Compute windows positions (hard coded for a 1920x1080 monitor)
window_pos1 = [4, 450, 400, 400];
window_pos2 = [1150, 450, 450, 400];
window_pos3 = [20, 50, 1100, 800];

%--- As a background independent process, plot the real-time EEG data and
% assign it to the windows position 1
% vis_stream();
% raw_data_figure_handle = gcf;
% set(raw_data_figure_handle, 'Position', window_pos1, ...
%     'MenuBar', 'none', 'ToolBar', 'none');

%% Fill sample buffer and pre-allocate figures

% Send event label (using LSL)
event_str = 'Buffer_start';
outlet.push_sample({event_str});

%--- Fill EEG buffer according to the input size (seconds)
[EEG_buffer, buffer_op_status] =  ...
    LSL_real_time_EEG_buffer_filling(sampling_rate, EEG_buffer_size_sec, number_of_chanels);

% Send event label (using LSL)
tic
event_str = 'Buffer_end';
outlet.push_sample({event_str});

%% --- Figure 2: EEG buffer Frequency Spectrum
% Use the Pre-processing function on the current EEG buffer
if buffer_op_status
    
    % Pre-process (filter, refer, etc.) EEG buffer
    [process_op_status, processed_EEG] = ...
        ssvepLiveProcessingAndDetection_v3(EEG_buffer, ...
        sampling_rate, band_pass_filter_freqs, common_reference);

    % Create freq. spectrum figure (invisible on the background)
    temp_fig = figure('Visible', 'off', 'Position', window_pos2);
    [freq_spectra_data, freq_values, ~, ~] = pop_spectopo(processed_EEG, 1, ...
        [processed_EEG.xmin (processed_EEG.xmax*1000)], ...
        'EEG' , 'freqrange', freq_range_spectra_vis, 'electrodes', 'off');
    delete(allchild(temp_fig));
    close(temp_fig);
    
    % Use spectrum data from figure and obtain values for freq. range
    [freq_pos1, freq_pos2] = getPositionsForValueRange(...
        freq_values, freq_range_spectra_vis);
    freq_values = freq_values(freq_pos1:freq_pos2);
    freq_data = freq_spectra_data(:, freq_pos1:freq_pos2);
    
    % Create the new (and visible) Frequency Spectrum plot
    freq_spectrum_figure_handle = figure('Position', window_pos2, ...
        'Renderer', 'OpenGL', 'MenuBar', 'none', 'ToolBar', 'none');
    
    % Add frequency spectrum plot to the figure
    freq_spectrum_axes_handle = axes('Parent', freq_spectrum_figure_handle);
    freq_spectrum_plot_handle = plot(freq_values, freq_data, 'LineWidth', 1);
    
    title(['EEG Buffer (', num2str(EEG_buffer_size_sec), ...
        ' Sec.) Frequency Spectrum'], 'FontSize', 18);
    
    ylabel('Power (units ommited)', 'FontSize', 16);
    grid on
    
    % legend(EEG_channel_10_20_labels);
    
    % Add dynamic xlabel with processing delay
    freq_spectrum_time = toc;
    freq_spectrum_xlabel_str = char(['Frequency (Hz).  [Processing Delay: ', ...
        num2str(freq_spectrum_time, 4), ' Sec.]']);
    
    freq_spectrum_xlabel_handle = xlabel(freq_spectrum_xlabel_str, ...
        'Parent', freq_spectrum_axes_handle, 'FontSize', 16);
    
    % Raise flag indicating the process ended correctly
    if process_op_status
        spectrum_op_status = true;
    else
        spectrum_op_status = false;
    end
end


%% --- Figure 4: frequency bar of each channel
tic;
% Perform CCA using the SSVEP stimulus frequencies
% [per_power] = freq_histgram(processed_EEG, freq_bank);

% Create (CCA's Rho values) figure, axes and plot
CCA_figure_handle = figure('Position', window_pos3, ...
    'Renderer', 'OpenGL', 'MenuBar', 'none', 'ToolBar', 'none');

xtick = 8:17;
freq_band = freq_values(xtick);
freq_data = freq_spectra_data(:, xtick);
xtick = xtick - xtick(1)+1;
xtick = xtick(1:2:end);
xtick_lbl = cell(1, length(xtick));
for idx = 1:length(xtick)
    xtick_lbl{idx} = num2str(round(freq_band(xtick(idx))));
end

bar_handles = zeros(1, channel_length);
for row = 1:8
    for col = 1:4
        axes_h = subplot(8,4,(row-1)*4+col, 'Parent', CCA_figure_handle);
        bar_h = bar(freq_data((row-1)*4+col, :),...
            'Parent', axes_h);
        bar_handles((row-1)*4+col) = bar_h;
        set(axes_h, ...
            'XTickLabel',xtick_lbl, ...
            'XTick',xtick);
        xlim(axes_h,[0 length(freq_band)+1]);
        ylim(axes_h, [-30 30]);
        title(processed_EEG.chanlocs((row-1)*4+col).labels);
    end
end

CCA_time = toc;
CCA_xlabel_str = char(['Stimulus Frequency (Hz).  [Processing Delay: ', ...
    num2str(freq_spectrum_time + CCA_time, 4), ' Sec.]']);
disp(CCA_xlabel_str);

%% Run simulation, updating figures everytime buffer is filled
while buffer_op_status
    
    % Send event label (using LSL)
    event_str = 'Buffer_start';
    outlet.push_sample({event_str});
    
    % Update EEG buffer
    [EEG_buffer, buffer_op_status] =  ...
        LSL_real_time_EEG_buffer_filling(sampling_rate, EEG_buffer_size_sec, number_of_chanels);
    
    % Send event label (using LSL)
    tic
    event_str = 'Buffer_end';
    outlet.push_sample({event_str});
    
    %% --- Update figure 2 (Frequency spectrum)
    if buffer_op_status
   
        % Pre-process (filter, refer, etc.) EEG buffer
        [process_op_status, processed_EEG] = ...
            ssvepLiveProcessingAndDetection_v3(EEG_buffer, ...
            sampling_rate, band_pass_filter_freqs, common_reference);
        
        % Create freq. spectrum figure (invisible on the background)
        temp_fig = figure('Visible', 'off', 'Position', window_pos2);
        [freq_spectra_data, freq_values, ~, ~] = pop_spectopo(processed_EEG, 1, ...
            [processed_EEG.xmin (processed_EEG.xmax*1000)], ...
            'EEG' , 'freqrange', freq_range_spectra_vis, 'electrodes', 'off');
        delete(allchild(temp_fig));
        close(temp_fig);
        
        % Update frequency spectrum info within visualization range
        freq_values = freq_values(freq_pos1:freq_pos2);
        freq_data = freq_spectra_data(:, freq_pos1:freq_pos2);
        
        % Update lines on Frequency Spectrum figure
        for id=1:length(freq_spectrum_plot_handle);
            set(freq_spectrum_plot_handle(id), 'Ydata', freq_data(id, :));
        end
        freq_spectrum_time = toc;
        
        % Update xlabel with processing delay
        freq_spectrum_xlabel_str = char(['Frequency (Hz).  [Processing Delay: ', ...
            num2str(freq_spectrum_time, 3), ' Sec.]']);
        
        set(freq_spectrum_xlabel_handle, 'String', freq_spectrum_xlabel_str);
        
        % Update figure
        drawnow expose
        
        % Raise flag indicating the process ended correctly
        if process_op_status
            spectrum_op_status = true;
        else
            spectrum_op_status = false;
        end
    end
    
    %% --- Update figure 4: CCA Score and frequency detection
        
    % Update figure's bar plot
    for row = 1:8
        for col = 1:4
            set(bar_handles((row-1)*4+col), 'Ydata', freq_data((row-1)*4+col, :));
        end
    end
    
    % Update Xlabel
    CCA_time = toc;
    CCA_xlabel_str = char(['Stimulus Frequency (Hz).  [Processing Delay: ', ...
        num2str(freq_spectrum_time + CCA_time, 4), ' Sec.]']);
    disp(CCA_xlabel_str);
    
    % Update figure
    drawnow expose
end