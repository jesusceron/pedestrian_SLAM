function [beac_rssi_fixed_filtered,beac_rssi_activity_filtered] = rssiKF(beac_rssi_fixed,beac_rssi_activity)

R_std = 10;
Q_std = 0.3;

%% beac_rssi_fixed
number_of_rows = size(beac_rssi_fixed, 1);
number_of_columns = size(beac_rssi_fixed, 2);

beac_rssi_fixed_filtered = zeros(number_of_rows, number_of_columns);

for i_column = 1:number_of_columns
    
    [sample_indexes, ~, rssis] = find(beac_rssi_fixed(:,i_column)); %find indexes with values different to zero
    rssis_filtered_py = py.mifunc.rssi_filter(rssis, R_std, Q_std); % Kalman Filter
    %class(rssis_filtered{1})
    rssis_filtered_cell = cell(rssis_filtered_py)'; % convert to double array
    rssis_filtered = cellfun(@double,rssis_filtered_cell); % save
    
    % To asign the rssis filtered values to the array
    for i_index = 1:size(sample_indexes)
        
        sample_index = sample_indexes(i_index);
        beac_rssi_fixed_filtered(sample_index, i_column) = rssis_filtered(i_index);
        
    end
    
end


%% beac_rssi_activity

number_of_rows = size(beac_rssi_activity, 1);
number_of_columns = size(beac_rssi_activity, 2);

beac_rssi_activity_filtered = zeros(number_of_rows, number_of_columns);

for i_column = 1:number_of_columns
    
    [sample_indexes, ~, rssis] = find(beac_rssi_activity(:,i_column)); %find indexes with values different to zero
    rssis_filtered_py = py.mifunc.rssi_filter(rssis, R_std, Q_std); % Kalman Filter
    %class(rssis_filtered{1})
    rssis_filtered_cell = cell(rssis_filtered_py)'; % convert to double array
    rssis_filtered = cellfun(@double,rssis_filtered_cell); % save
    
    % To asign the rssis filtered values to the array
    for i_index = 1:size(sample_indexes)
        
        sample_index = sample_indexes(i_index);
        beac_rssi_activity_filtered(sample_index, i_column) = rssis_filtered(i_index);
        
    end
    
end


end
