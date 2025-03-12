function [corr_left_neuron, corr_right_neuron] = process_trials_neural_corr(chosen_mice, current_dataset, mouse_date, dff_st, frames_before, frames_after, trials_left, trials_right,combined_thres,current_cell_type, function_params, save_data_directory)
        
    % LOAD ALIGNED NEURAL DATA
    field_name_neural = function_params.field_name_neural;
    aligned_neural_data = dff_st{1, current_dataset}.(field_name_neural);
    data_norm = normalize_aligned_data(aligned_neural_data, 'zscore', []);

    %get specified cell type
    if isempty(combined_thres)
        cellCount = size(dff_st{1,current_dataset}.(field_name_neural),2); %loading from control which does not have opto so only sounds!
        mod_cells = 1:cellCount;
        sig_mod_boot{1,current_dataset} = mod_cells;
    else
        sig_mod_boot{1,current_dataset} = combined_thres{1,current_dataset};
    end
    
    left_matrix = data_norm(trials_left, sig_mod_boot{1,current_dataset}, [frames_before,frames_after]);
    right_matrix = data_norm(trials_right, sig_mod_boot{1,current_dataset}, [frames_before,frames_after]);

    %comparison matrix
    comparison_cells = setdiff(function_params.current_celltypes.(current_cell_type),sig_mod_boot{1,current_dataset}); %do not correlate with itself
    left_matrix_comp = data_norm(trials_left, comparison_cells, [frames_before,frames_after]);
    right_matrix_comp = data_norm(trials_right, comparison_cells, [frames_before,frames_after]);

    

% Binning option
    if function_params.corr_bins > 0
        bin_size = function_params.corr_bins;
    
        if function_params.nonoverlap_bins == 1
            % Non-overlapping binning
            num_bins = floor(size(left_matrix, 3) / bin_size);
                        
            % Reshape and take mean for neural data
            left_matrix = squeeze(mean(reshape(left_matrix(:, :, 1:num_bins * bin_size), [], size(left_matrix, 2), bin_size, num_bins), 3));
            right_matrix = squeeze(mean(reshape(right_matrix(:, :, 1:num_bins * bin_size), [], size(right_matrix, 2), bin_size, num_bins), 3));

            left_matrix_comp = squeeze(mean(reshape(left_matrix_comp(:, :, 1:num_bins * bin_size), [], size(left_matrix_comp, 2), bin_size, num_bins), 3));
            right_matrix_comp = squeeze(mean(reshape(right_matrix_comp(:, :, 1:num_bins * bin_size), [], size(right_matrix_comp, 2), bin_size, num_bins), 3));

            
        else
            % Overlapping binning using movmean            
            left_matrix = movmean(left_matrix, bin_size, 3, 'Endpoints', 'discard');
            right_matrix = movmean(right_matrix, bin_size, 3, 'Endpoints', 'discard');

            left_matrix_comp = movmean(left_matrix_comp, bin_size, 3, 'Endpoints', 'discard');
            right_matrix_comp = movmean(right_matrix_comp, bin_size, 3, 'Endpoints', 'discard');

        end
    end

    
    % CONCATENATE TRIALS HORIZONTALLY    
    left_concatenated_neural = reshape(permute(left_matrix, [3, 1, 2]), [], size(left_matrix, 2));
    right_concatenated_neural = reshape(permute(right_matrix, [3, 1, 2]), [], size(right_matrix, 2));

    left_concatenated_neural_comp = reshape(permute(left_matrix_comp, [3, 1, 2]), [], size(left_matrix_comp, 2));
    right_concatenated_neural_comp = reshape(permute(right_matrix_comp, [3, 1, 2]), [], size(right_matrix_comp, 2));

    
    for neuron = 1:size(left_concatenated_neural, 2)
        corr_vals_left = zeros(1, size(left_concatenated_neural_comp, 2));
        corr_vals_right = zeros(1, size(right_concatenated_neural_comp, 2));
        
        for comp_neuron = 1:size(left_concatenated_neural_comp, 2)
            [r, p] = corrcoef(left_concatenated_neural(:, neuron), left_concatenated_neural_comp(:, comp_neuron));
            corr_vals_left(comp_neuron) = r(1, 2);
            
            [r, p] = corrcoef(right_concatenated_neural(:, neuron), right_concatenated_neural_comp(:, comp_neuron));
            corr_vals_right(comp_neuron) = r(1, 2);
        end
        
        % Take mean correlation across all comparison neurons
        corr_left_neuron(neuron) = mean(corr_vals_left, 'omitnan');
        corr_right_neuron(neuron) = mean(corr_vals_right, 'omitnan');
    end

    
    % SAVE DATA IF NEEDED
    if ~isempty(save_data_directory)
        save(fullfile(save_data_directory, sprintf('%s_%s_%s_neural_corr.mat', chosen_mice{current_dataset}, mouse_date{current_dataset}, current_cell_type)), 'corr_left_neuron', 'corr_right_neuron', 'p_values_left', 'p_values_right');
    end
end
