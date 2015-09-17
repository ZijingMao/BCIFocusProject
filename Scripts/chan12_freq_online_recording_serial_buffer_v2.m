% Version 4, updated help will be added later
clear
close all
clc
zz=0;
aaa=1;
record_att_flag=1;
buf_att=1;
load('threshold.mat');%%the value 

record_att=0;
record_med=0;

% Add EEGlab and needed file\folders to the path
% exe_folder = 'E:\ProgramData\Dropbox\EEG Projects\Realtime\';  % Specific recording folder to use
exe_folder = 'G:\Dropbox\EEG Projects\Realtime\';  % Specific recording folder to use

cd([exe_folder 'Scripts']);
use_parfor = false;                      % active or not parallel for loops
[dataset_folder, script_folder] = set_path_and_parpool(use_parfor, exe_folder);

%% SSVEP frequency detection fine-tunning
CCA_num_harmonics = 2;
% EEG_buffer_size_sec = 3;          % EEG buffer Length (adjust per participant)
EEG_buffer_size_sec = 1;          % EEG buffer Length (adjust per participant)
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
channel_length = 12;
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
sampling_rate = 540;              % Fron BioSemi device OR simulation
number_of_chanels = length(EEG_channel_10_20_labels); % Update this to match experiment

%--- Compute windows positions (hard coded for a 1920x1080 monitor)
window_pos1 = [4, 450, 400, 400];
window_pos2 = [1150, 450, 450, 400];
window_pos3 = [20, 50, 1100, 800];
window_pos4 = [800, 450, 500, 500];

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
        sampling_rate, band_pass_filter_freqs, common_reference, 12, 1);

    % Create freq. spectrum figure (invisible on the background)
    temp_fig = figure('Visible', 'off', 'Position', window_pos2);
    [freq_spectra_data, freq_values, ~, ~] = pop_spectopo(processed_EEG, 1, ...
        [processed_EEG.xmin (processed_EEG.xmax*1000)], ...
        'EEG' , 'freqrange', freq_range_spectra_vis, 'electrodes', 'off');%%Here makes the figure 1 tempfig
    delete(allchild(temp_fig));
    close(temp_fig);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
for row = 1:3
    for col = 1:4
        axes_h = subplot(3,4,(row-1)*4+col, 'Parent', CCA_figure_handle);
        bar_h = bar(freq_data((row-1)*4+col, :),...
            'Parent', axes_h);
        bar_handles((row-1)*4+col) = bar_h;
        set(axes_h, ...
            'XTickLabel',xtick_lbl, ...
            'XTick',xtick);
        xlim(axes_h,[0 length(freq_band)+1]);
        ylim(axes_h, [-30 30]);
        if ~isempty(processed_EEG.chanlocs)
            title(processed_EEG.chanlocs((row-1)*4+col).labels);
        end
    end
end
%%%%%%%%%%%%%%%here comes figure 2

CCA_time = toc;
CCA_xlabel_str = char(['Stimulus Frequency (Hz).  [Processing Delay: ', ...
    num2str(freq_spectrum_time + CCA_time, 4), ' Sec.]']);
disp(CCA_xlabel_str);
%% attention and meditation level bar 
m=[5 0; 0 5];
  bar_handle = createfigure(m);
   ylim([1 10]) 
    title('attention and meditation level','FontSize', 18);
    
    ylabel('level', 'FontSize', 16);


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
            sampling_rate, band_pass_filter_freqs, common_reference, 12, 1);
        
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
   %%       %%%%%%%try to collect EEG data from realtime
    [mm nn]=size(freq_data);
    
    for iii=1:mm
        ff(zz*mm+iii,:)=freq_data(iii,:);
    end
    
    zz=zz+1;

%%
delta=freq_data(:,1:3);
theta=freq_data(:,4:7);
alpha=freq_data(:,8:12);
beta=freq_data(:,13:23);
aa=0;bb=0;cc=0;dd=0; p_delta=0; p_beta=0;p_theta=0;p_ratio_thetaoveralpha=0; p_ratio_thetaalphaoverbeta=0;
for i=1:12 %%12 channel
  
%    ratio(i)=ave_theta(i)/ave_alpha(i);%%%ratio
p_delta(i)=10^(delta(i)/10);
p_beta(i)=10^(beta(i)/10);
p_theta(i)=10^(theta(i)/10);
    p_alpha(i)=10^(alpha(i)/10);
    p_ratio_thetaoveralpha(i)=p_theta(i)/p_alpha(i);
%     p_ratio_betaoveralpha(i)=p_beta(i)/p_alpha(i);
    p_ratio_thetaalphaoverbeta(i)=(p_theta(i)+p_alpha(i))/p_beta(i);
end
for k=1:4
    alpha_rate(k)=(alpha(k)-threshold_alpha(2,k))/(threshold_alpha(1,k)-threshold_alpha(2,k));
    beta_rate(k)=(beta(k)-threshold_beta(2,k))/(threshold_beta(1,k)-threshold_beta(2,k));
    theta_rate(k)=(theta(k)-threshold_theta(2,k))/(threshold_theta(1,k)-threshold_theta(2,k));   
    ratio_thetaalphaoverbeta(k)=(p_ratio_thetaalphaoverbeta(k)-ratio_thetaalphabeta(2,k))/(ratio_thetaalphabeta(1,k)-ratio_thetaalphabeta(2,k));
    ratio_thetaoveralpha(k)=(p_ratio_thetaoveralpha(k)-ratio_thetaalpha(2,k))/(ratio_thetaalpha(1,k)-ratio_thetaalpha(2,k));
end

for m=6:9
    ratio_thetaalphaoverbeta(m-1)=(p_ratio_thetaalphaoverbeta(m)-ratio_thetaalphabeta(2,m-1))/(ratio_thetaalphabeta(1,m-1)-ratio_thetaalphabeta(2,m-1));
    ratio_thetaoveralpha(m-1)=(p_ratio_thetaoveralpha(m)-ratio_thetaalpha(2,m-1))/(ratio_thetaalpha(1,m-1)-ratio_thetaalpha(2,m-1));

end
for j=10:12

 alpha_rate(j-5)=(alpha(j)-threshold_alpha(2,j-5))/(threshold_alpha(1,j-5)-threshold_alpha(2,j-5));
    beta_rate(j-5)=(beta(j)-threshold_beta(2,j-5))/(threshold_beta(1,j-5)-threshold_beta(2,j-5));
    theta_rate(j-5)=(theta(j)-threshold_theta(2,j-5))/(threshold_theta(1,j-5)-threshold_theta(2,j-5));
     ratio_thetaalphaoverbeta(j-1)=(p_ratio_thetaalphaoverbeta(j)-ratio_thetaalphabeta(2,j-1))/(ratio_thetaalphabeta(1,j-1)-ratio_thetaalphabeta(2,j-1));
    ratio_thetaoveralpha(j-1)=(p_ratio_thetaoveralpha(j)-ratio_thetaalpha(2,j-1))/(ratio_thetaalpha(1,j-1)-ratio_thetaalpha(2,j-1));

end
for k=1:7
if theta_rate(k)>1 
    theta_rate(k)=1;
end
    if theta_rate(k)<0
        theta_rate(k)=0;
    end
    
if alpha_rate(k)>1 
    alpha_rate(k)=1;
end
    if alpha_rate(k)<0
        alpha_rate(k)=0;
    end
    
    if  beta_rate(k)>1 
     beta_rate(k)=1;
end
    if  beta_rate(k)<0
         beta_rate(k)=0;
    end
    med(k)=beta_rate(k)/3+alpha_rate(k)/3+theta_rate(k)/3;
end

for k=1:11
    if ratio_thetaalphaoverbeta(k)>1 
    ratio_thetaalphaoverbeta(k)=1;
    end
    if ratio_thetaalphaoverbeta(k)<0
        ratio_thetaalphaoverbeta(k)=0;
    end
    if ratio_thetaoveralpha(k)>1 
    ratio_thetaoveralpha(k)=1;
    end
    if ratio_thetaoveralpha(k)<0
        ratio_thetaoveralpha(k)=0;
    end
    att(k)=ratio_thetaalphaoverbeta(k)*0.5+ ratio_thetaoveralpha(k)*0.5;
end
att_lvl=mean(att);
med_lvl=mean(med);

aaa=aaa+1;

record_att_1=[record_att round(att_lvl*10)];
record_med_1=[record_med round(med_lvl*10)];

buf_att=[buf_att round(att_lvl*10)];

if att_lvl>0.5
    att_flag=1;
    if att_lvl<=0.5
        att_flag=0;
    end
end

record_att_flag=[record_att_flag att_flag];


if record_att_flag(aaa)==1 && record_att_flag(aaa-1)==1
    buf_att(aaa)=buf_att(aaa-1);
   
end 

if record_att_flag(aaa)==0 && record_att_flag(aaa-1)==1
    buf_att=record_att_flag(aaa-1);
end 

if record_att_flag(aaa)==0 && record_att_flag(aaa-1)==0
    buf_att=0;
end 

if record_att_flag(aaa)==1 && record_att_flag(aaa-1)==0
    buf_att=record_att_flag(aaa-1);
end 


    

%%

        
        % Update lines on Frequency Spectrum figure
        for id=1:length(freq_spectrum_plot_handle);
            set(freq_spectrum_plot_handle(id), 'Ydata', freq_data(id, :));
%             get(freq_spectrum_plot_handle(id),'Ydata')
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
    for row = 1:3
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
    %% updata bar 
    for k=1:length(bar_handle)
%         mm=[rand(1)*10 0; 0 rand(1)*10];
        mm=[round(att_lvl*10) 0; 0 round(med_lvl*10)];
        set(bar_handle(k),'Ydata', mm(k,:));

        end
        drawnow expose
    
end