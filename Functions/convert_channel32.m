function [EEG] = convert_channel32(EEG)

chanlocs = EEG.chanlocs;

[~,~,Electrode] = xlsread('cog_system_32.xlsx', 'A1:A32');
M = xlsread('cog_system_32.xlsx');
for idx = 1:32
    for idx1 = 1:32
        currentchan = chanlocs(idx);
        if strcmpi(currentchan.labels, Electrode{idx1})
            currentchan.theta = M(idx1, 1);
            currentchan.radius = M(idx1, 2);
            currentchan.X = -M(idx1, 3);
            currentchan.Y = -M(idx1, 4);
            currentchan.Z = M(idx1, 5);
        end
        chanlocs(idx) = currentchan;
    end
end

EEG.chanlocs = chanlocs;

end
