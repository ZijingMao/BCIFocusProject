function [per_power] =  freq_histgram(EEG_data, ssvep_freq)

% Obtain information from the EEG data
sampling_rate = EEG_data.srate;

% Generate SSVEP stimulus signals]
per_power = zeros(EEG_data.nbchan, length(ssvep_freq)-1);
for id=1:length(ssvep_freq)-1
    
    freq_band = ssvep_freq(id):ssvep_freq(id+1);
    for idx_ch = 1:EEG_data.nbchan
        pband = bandpower(EEG_data.data(idx_ch,:),sampling_rate,freq_band);
        per_power(idx_ch, id) = pband;
    end
    
end

end
