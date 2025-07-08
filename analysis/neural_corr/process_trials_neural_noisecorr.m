function [corr_left_neuron, corr_right_neuron] = process_trials_neural_noisecorr(chosen_mice, current_dataset, mouse_date, dff_st, frames_before, frames_after, trials_left, trials_right,combined_thres,current_cell_type, function_params, save_data_directory)
        
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
        mod_cells = combined_thres{1,current_dataset};
    end
    comparison_cells = setdiff(function_params.current_celltypes.(current_cell_type), mod_cells); % don't correlate with self

    
    left_data = data_norm(trials_left, sig_mod_boot{1,current_dataset}, [frames_before,frames_after]);
    right_data = data_norm(trials_right, sig_mod_boot{1,current_dataset}, [frames_before,frames_after]);

    % === SUBTRACT MEAN ACROSS TRIALS TO GET NOISE ===
    mean_left = squeeze(mean(left_data, 1));   % [neurons x time]
    mean_right = squeeze(mean(right_data, 1)); % [neurons x time]

    noise_left = bsxfun(@minus, left_data, permute(mean_left, [3, 1, 2]));    % [trials x neurons x time]
    noise_right = bsxfun(@minus, right_data, permute(mean_right, [3, 1, 2]));

    % === BINNING ===
    if function_params.corr_bins > 0
        bin_size = function_params.corr_bins;

        if function_params.nonoverlap_bins == 1
            num_bins = floor(size(noise_left, 3) / bin_size);
            noise_left = squeeze(mean(reshape(noise_left(:, :, 1:num_bins * bin_size), [], size(noise_left, 2), bin_size, num_bins), 3));
            noise_right = squeeze(mean(reshape(noise_right(:, :, 1:num_bins * bin_size), [], size(noise_right, 2), bin_size, num_bins), 3));
        else
            noise_left = movmean(noise_left, bin_size, 3, 'Endpoints', 'discard');
            noise_right = movmean(noise_right, bin_size, 3, 'Endpoints', 'discard');
        end
    else
        noise_left = squeeze(mean(noise_left, 3));   % average across time
        noise_right = squeeze(mean(noise_right, 3));
    end

    % === RESHAPE TO [samples x neurons] ===
    left_matrix = reshape(permute(noise_left, [3 1 2]), [], size(noise_left, 2));   % [samples x neurons]
    right_matrix = reshape(permute(noise_right, [3 1 2]), [], size(noise_right, 2));

    % Separate modulated and comparison neurons
    left_mod = left_matrix(:, mod_cells);
    left_comp = left_matrix(:, comparison_cells);

    right_mod = right_matrix(:, mod_cells);
    right_comp = right_matrix(:, comparison_cells);

    % === COMPUTE NOISE CORRELATIONS ===
    corr_left_neuron = nan(1, size(left_mod, 2));
    corr_right_neuron = nan(1, size(right_mod, 2));

    for i = 1:size(left_mod, 2)
        corrs = corr(left_mod(:, i), left_comp, 'rows', 'pairwise');
        corr_left_neuron(i) = mean(corrs, 'omitnan');
    end

    for i = 1:size(right_mod, 2)
        corrs = corr(right_mod(:, i), right_comp, 'rows', 'pairwise');
        corr_right_neuron(i) = mean(corrs, 'omitnan');
    end

    % === SAVE IF NEEDED ===
    if ~isempty(save_data_directory)
        save(fullfile(save_data_directory, sprintf('%s_%s_%s_noise_corr.mat', ...
            chosen_mice{current_dataset}, mouse_date{current_dataset}, current_cell_type)), ...
            'corr_left_neuron', 'corr_right_neuron');
    end
end



