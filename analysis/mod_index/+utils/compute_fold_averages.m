function [data1_fold_avg, data2_fold_avg] = compute_fold_averages(data_fold,data2_fold, response_range, mod_type)
%compute average across frames for data!

% Compute averages based on mod_type
if strcmp(mod_type, 'prepost') && length(response_range) > 1 
    data1_fold_avg = mean(data_fold(:, :, response_range{1}), 3); %use first only (assuming stim)
    data2_fold_avg = mean(data_fold(:, :, response_range{2}), 3); %use first only (assuming stim)
elseif strcmp(mod_type, 'prepost_sound') && length(response_range) > 1
    data1_fold_avg = mean(data2_fold(:, :, response_range{1}), 3); %use first only (assuming ctrl) for sounds
    data2_fold_avg = mean(data2_fold(:, :, response_range{2}), 3); %use first only (assuming ctrl) for sounds
elseif strcmp(mod_type, 'prepost_ctrl')
    data1_fold_avg = mean(data_fold(:, :, response_range{1}), 3) - mean(data_fold(:, :, response_range{2}), 3); %difference of values pre-post
    data2_fold_avg = mean(data2_fold(:, :, response_range{1}), 3)- mean(data2_fold(:, :, response_range{2}), 3);
else
    data1_fold_avg = mean(data_fold(:, :, response_range{1}), 3);
    data2_fold_avg = mean(data2_fold(:, :, response_range{1}), 3);
end
end