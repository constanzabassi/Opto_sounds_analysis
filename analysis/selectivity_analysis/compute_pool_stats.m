function stats = compute_pool_stats(pool_data)
    % Compute summary statistics for a pool
    stats = struct();
    stats.mean_left = mean(pool_data.left_mod);
    stats.mean_right = mean(pool_data.right_mod);
    stats.std_left = std(pool_data.left_mod);
    stats.std_right = std(pool_data.right_mod);
    stats.num_cells = length(pool_data.left_mod);
    stats.percent_left_pref = sum(strcmp(pool_data.preferred_side, 'left')) / stats.num_cells * 100;
end