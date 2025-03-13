function [snr_left_avg, snr_right_avg, snr_left, snr_right] = process_trials_neural_snr(chosen_mice, current_dataset, mouse_date, dff_st, frames_before, frames_after, trials_left, trials_right, combined_thres,  function_params, save_data_directory)

    % LOAD ALIGNED NEURAL DATA
    field_name_neural = function_params.field_name_neural;
    aligned_neural_data = dff_st{1, current_dataset}.(field_name_neural);
    data_norm = normalize_aligned_data(aligned_neural_data, 'zscore', []);

    % Get specified cell type
    if isempty(combined_thres)
        cellCount = size(dff_st{1, current_dataset}.(field_name_neural), 2);
        mod_cells = 1:cellCount;
        sig_mod_boot{1, current_dataset} = mod_cells;
    else
        sig_mod_boot{1, current_dataset} = combined_thres{1, current_dataset};
    end

    % Extract neural data for significant neurons
    left_matrix = data_norm(trials_left, sig_mod_boot{1, current_dataset}, [frames_before, frames_after]);
    right_matrix = data_norm(trials_right, sig_mod_boot{1, current_dataset}, [frames_before, frames_after]);

    % Binning option
    if function_params.corr_bins > 0
        bin_size = function_params.corr_bins;

        if function_params.nonoverlap_bins == 1
            num_bins = floor(size(left_matrix, 3) / bin_size);

            left_matrix = squeeze(mean(reshape(left_matrix(:, :, 1:num_bins * bin_size), [], size(left_matrix, 2), bin_size, num_bins), 3));
            right_matrix = squeeze(mean(reshape(right_matrix(:, :, 1:num_bins * bin_size), [], size(right_matrix, 2), bin_size, num_bins), 3));
        else
            left_matrix = movmean(left_matrix, bin_size, 3, 'Endpoints', 'discard');
            right_matrix = movmean(right_matrix, bin_size, 3, 'Endpoints', 'discard');
        end
    end


    % Compute SNR for significant neurons (Signal-to-Noise Ratio)

    mean_response_left = squeeze(mean(left_matrix, 1));  % Average across trials (dim 1)
    std_response_left = squeeze(std(left_matrix, 0, 1)); % STD across trials (dim 1)
    snr_left = mean_response_left ./ std_response_left; %neurons x time bins
    snr_left_avg = mean(snr_left, 2); % Averaging across time bins (dim 2)

    mean_response_right = squeeze(mean(right_matrix, 1));  % Average across trials (dim 1)
    std_response_right = squeeze(std(right_matrix, 0, 1)); % STD across trials (dim 1)
    snr_right = mean_response_right ./ std_response_right;
    snr_right_avg = mean(snr_right, 2); % Averaging across time bins (dim 2)


    % SAVE DATA IF NEEDED
%     if ~isempty(save_data_directory)
%         save(fullfile(save_data_directory, sprintf('%s_%s_neural_snr.mat', chosen_mice(current_dataset), mouse_date{current_dataset})), ...
%             'snr_left', 'snr_right');
%     end
end
