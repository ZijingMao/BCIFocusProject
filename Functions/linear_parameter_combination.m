function parameter_combination = ...
    linear_parameter_combination(input_parameteres_cell_array)
%LINEAR_PARAMETER_COMBINATION creates an array of every combination of
%input parameters, each set of input as a cell on the input cell array.
%Each parameter should contain a vector with the set of values used on each
%variable.

% Pre-allocate output
parameter_combination = cell(1, length(input_parameteres_cell_array));

% Get the number of parameter values to combine
parameters_set_length = cellfun(@length, input_parameteres_cell_array);
parameters_set_size_product = prod(parameters_set_length);

% Fill each output array
for param_num=1:length(input_parameteres_cell_array)
    
    % Pre-allocate array on the output cell array
    %-- Numeric parameter ser
    if sum(isdouble(input_parameteres_cell_array{param_num})) ...
            == length(input_parameteres_cell_array{param_num})
        
        parameter_combination{param_num} = zeros(1, parameters_set_size_product);
        numeric_flag = true;
    else
        
        % String, structs, cells, mixed. Create an inner cell array
        parameter_combination{param_num} = cell(1, parameters_set_size_product);
        cell_flag = true;
    end
    
    
    % Fill the current parameter set
    element_switch_counter = 2^(param_num - 1);
    current_parameter_values = input_parameteres_cell_array{param_num};
    
    % Loop through    
    
    
end


end

