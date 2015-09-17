function cca_corr_rho = ssvepCCAcorr(EEG_data, SSVEP_freq)
% This functions computes the canonical correlation of some input
% multi-channel EEG data (from SSVEP experiment mainly but any kind of EEG
% experiment can be entered as input) and the first harmonic of some input
% frequencies [sin(2*pi*f*N)] where N is the number of samples of the input
% EEG data, these are the SSVEP stimulus frequencies. Input EEG data should
% be (at least for now) on current EEGLab format. On CCA algorithm by
% default 6 harmonics are used.

% Obtain information from the EEG data
time_points = 1:length(EEG_data.data);
sampling_rate = EEG_data.srate;

% Generate SSVEP stimulus signals
num_harmonics = 5;
cca_Y = cell(1, length(SSVEP_freq));
for id=1:length(SSVEP_freq)
    
    % Pre-allocate Yf for current freq.
    cca_Y{id} = zeros(num_harmonics*2, length(time_points));
    
    % Generate reference signal matrices Y's
    N_ = time_points ./ sampling_rate;
    for id1=1:num_harmonics
        for id2=1:length(EEG_data.data)
            cca_Y{id}( (2*(id1-1)+1), id2) = sin( 2*pi*id1*SSVEP_freq(id)*N_(id2) );
            cca_Y{id}( 2*id1, id2)         = cos( 2*pi*id1*SSVEP_freq(id)*N_(id2) );
        end
    end
end

% Test CCA for prediction of the input data
cca_corr_rho = zeros(1, length(SSVEP_freq));
for id=1:length(SSVEP_freq)
    [~, ~, rho_, ~, ~, ~] = canoncorr(double(EEG_data.data)', cca_Y{id}');
    cca_corr_rho(id) = rho_(1);
end

% Plot result as an intuitive graph IF option selected
% if plot_flag == 1
%     figure, bar([R1;R2], 'grouped'); grid on
%     title(['Canonical Correlation (r), 2 samples, ', num2str(num_channels), ...
%         ' Channels'], 'Fontsize', 16);
%     xlabel('Simulated Samples', 'Fontsize', 14);
%     ylabel('Correlation (the closer to 1 the better?)', 'Fontsize', 14);
%     hx = legend(chan_str);
%     set(hx, 'FontSize', 12);
%     
%     figure, subplot(1, 2, 1);
%     plot(V1, U1, '.', 'LineWidth', 2); grid on
%     axis([min(reshape([V1;V2], 1, [])) max(reshape([V1;V2], 1, []))...
%         min(reshape([U1;U2], 1, [])) max(reshape([U1;U2], 1, []))]);
%     title('Sample X1', 'Fontsize', 16);
%     xlabel('Canonical Scores V', 'FontSize', 14);
%     ylabel('Canonical Scores U', 'FontSize', 14);
%     legend(chan_str, 'Location', 'Best');
%     
%     subplot(1, 2, 2);
%     plot(V2, U2, '.', 'LineWidth', 2); grid on
%     axis([min(reshape([V1;V2], 1, [])) max(reshape([V1;V2], 1, []))...
%         min(reshape([U1;U2], 1, [])) max(reshape([U1;U2], 1, []))]);
%     title('Sample X2', 'Fontsize', 16);
%     xlabel('Canonical Scores V', 'FontSize', 14);
%     ylabel('Canonical Scores U', 'FontSize', 14);
%     legend(chan_str, 'Location', 'Best');
% end


end
