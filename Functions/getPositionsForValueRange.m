function [start_pos, end_pos] = getPositionsForValueRange(...
    array_values, array_target_range)
% Simpleistic function to get the array positions that correspond to a
% range of values (of course present on the array)

start_pos = (1:length(array_values)) .* (array_values >= array_target_range(1))';
start_pos( start_pos == 0 ) = [];
start_pos = min(start_pos);

end_pos = (1:length(array_values)) .* (array_values <= array_target_range(2))';
end_pos = max(end_pos);


end

