function [op_status, EEG] = ssvepLiveProcessingAndDetection_v3(...
    EEG_epoch, Fs, pb_filter, common_reference, channel, chanOffset)

% This function receives a chunk of EEG data from an SSVEP experiment,
% recorded using either OFFLINE or ONLINE mode. The function performs
% essential pre-processing (EEGLAB-based, already in path) such:
%   - Downsampling to 512 samples/second
%   - FIR filtering (pass-band, cut-off freqs as inputs)
%   - Re-reference by channel average
%   - Baseline removal (data is a single epoch)
%   - Data frequency spectra (optional, plotted on a window or graphic data
%     extracted as an extra output)
%
% When the input data epoch is pre-processed, the the SSVEP frequency
% detection is done by Canonical Correlation Analysis (CCA) and based on
% its rho output a selection is from the input SSVEP frequencies.
%
%  INPUTS:
%   - EEG epoch: an OFFLINE or ONLINE recorded piece of data (or obtained in
%     real-time from LSL library).
%   - Fs: Sampling rate of the input data.
%   - pb_filter: a 1x2 double array with the cut-off frequencies for the
%     pass-band FIR filter, if the SSVEP presentation program is displayed
%     on a standard 60Hz monitor, the recommended values are [1, 60].
%   - SSVEP_freqs: this are the frequencies for the SSVEP flickering
%     visual stimulus as displayed on the monitor, the values are required
%     to accurately determine which stimulus was observed based on the CCA
%     output.
%
%  OUTPUTS:
%   - op_status: a simple flag to indicated that the whole process was
%     completed sucessfully or if there was an input error, execution error,
%     etc. Values: {'OK', 'Input Error', 'Execution Error'}
%   - freq_decision: this is one the input SSVEP frequency values, the one
%     that according to the classifier is the most like to be present on the
%     input EEG data.
%   - rho_values: the actual output of the CCA prediction.
%   - EEG: an EEGLAB-like structure with the information obtained from
%     input and some other computed during the function execution.
%
%  USAGE:
%   - predicted frequency on 10Hz and 15Hz stimulus data, filtering 1-60Hz
%    [op_status, freq_decision, ~, ~] =
%    ssvepLiveProcessingAndDetection_v1(EEG_epoch, [1, 60], [10, 15], 'none')
%
%   - General form:
%    [op_status, freq_decision, rho_values, freq_spectrum_data] = ...
%    ssvepLiveProcessingAndDetection_v1(EEG_epoch, pb_filter, SSVEP_freqs, spectrum_flag)


%--- Build and EEGLAB-like input struct.
EEG.srate = Fs;
EEG.xmin = 0;
baseline_adjust_samples = 5;
EEG_epoch(end-chanOffset:end, :) = [];
EEG.xmax = (size(EEG_epoch, 2) - baseline_adjust_samples)/Fs;
EEG.trials = 1;
EEG.setname = 'xdf file';
EEG.chaninfo = struct('plotrad', {[]}, 'shrink', {[]}, ...
    'nosedir', {'+X'}, 'nodatchans', {[]}, 'icachansind', {[]});

% No ICA performed, so fields are added empty
EEG.event = [];
EEG.icawinv = [];
EEG.icaweights = [];
EEG.icasphere = [];
EEG.icaact = [];
EEG.icachansind = [];

%--- Perform Re-referencing based on all channels
% ref_mean = repmat(mean(EEG_epoch), [size(EEG_epoch, 1), 1]);
% EEG_epoch = EEG_epoch - ref_mean;

%--- Remove bad channels, extract recorded data from epoch
%EEG_epoch = EEG_epoch(EEG_LSL_matrix_position, :);
EEG_epoch( isnan(EEG_epoch) ) = 0;

%--- Perform pre-processing on input EEG data. It is assumed that input
% data is previously cleaned from empty and NaN values.
% Use generic chanlocs file that only contains labels A1 = Axx <= 32
% this is a temporal workaround for SSVEP recordings only, include proper
% and accurate channel locations as a required input soon.

% Load 32-channel sample locations
chanlocs = [];
fileName = ['EEG_channel_locations_' num2str(channel) '_cog_system.mat'];
if exist(fileName, 'file')
    load();
    EEG.chanlocs = chanlocs;
else
    EEG.chanlocs = [];
end
%EEG.chanlocs = sample_channels(EEG_channel_10_20_position);
EEG.data = EEG_epoch;
EEG.pnts = size(EEG.data, 2);
EEG.nbchan = size(EEG.data, 1);

try
    
    % Pass-band FIR filtering, Re-referencing (by channel average)
    EEG = pop_eegfiltnew(EEG, [], pb_filter(1), [], true, [], 0);
    EEG = pop_eegfiltnew(EEG, [], pb_filter(2), [], 0, [], 0);
    
    % Apply Common re-referencing if selected (At least have 4 channels)
    if common_reference && (EEG.nbchan > 3)
        EEG = pop_reref(EEG, []);
        %EEG.data = double(EEG.data);  % Possible fix to OFFLINE Problem???
    end
    
    % Baseline removal and remove reference channels data
    EEG = pop_rmbase(EEG, [EEG.xmin (EEG.xmax*1000)]);
    op_status = true;
    
catch
    
    % Here, something were wrong during the pre-processing,
    % report to the user and exit.
    disp('An ERROR has been produced during the process execution, please verify inputs');
    op_status = false;
    return;
end

end