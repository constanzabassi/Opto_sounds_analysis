function mod_index_stats = plot_context_comparisons(contexts_to_compare,overlap_labels, mod_indexm, sig_mod_boot, sorted_cells, params, save_dir)
    % Plot all context comparison visualizations
    
    % Create output directory if it doesn't exist
%     save_dir = params.mod.savepath;
    if ~exist(save_dir, 'dir') && ~isempty(save_dir)
        mkdir(save_dir);
    end


    % 2. Plot overlap between contexts
%     contexts_to_compare = [1,2]; % Active vs Passive
%     overlap_labels = {'Active', 'Passive', 'Both'};
    percent_cells = calculate_sig_overlap_percentages(...
        sig_mod_boot, ...
        mod_indexm, contexts_to_compare);
    plot_sig_overlap_pie(percent_cells, overlap_labels, save_dir, contexts_to_compare);

    % 3. Plot modulation index distributions
    [context_mod_all, celltypes_ids] = organize_sig_mod_index_contexts_celltypes(...
        params, mod_indexm, sig_mod_boot, sorted_cells);
    
    % Generate visualization suite
    [mod_index_stats] = generate_mod_index_plots(context_mod_all, celltypes_ids, params, save_dir);
end

