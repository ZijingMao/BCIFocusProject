function [ per_power ] = bandpower_cal( x, freq_band, sampling_rate )
%EEGSOLUTIONS_BANDPOWER Summary of this function goes here
%   Detailed explanation goes here

Freq = sampling_rate;
[~, channels, epoch] = size(x);

per_power = zeros(channels, epoch);

for idx = 1 : epoch

    for idx_ch = 1:channels
        pband = bandpower(x(:, idx_ch, idx),Freq,freq_band);
        per_power(idx_ch, idx) = pband;
    end

end

end

