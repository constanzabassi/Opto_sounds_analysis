function wrapper_selecitivity_pool_analysis(base, params, mod_indexm, sig_mod_boot, mod_index_results, avg_results, sorted_cells, all_celltypes, selectivity_indexm, data_type, varargin)
    % WRAPPER FUNCTION FOR SELECTIVITY POOL ANALYSIS
    %
    % Inputs:
    %   base - Base directory for saving results
    %   params - Structure containing analysis parameters
    %   mod_indexm - Modulation index matrix
    %   sig_mod_boot -Indices of modulated neurons (datasets,contexts)
    %   mod_index_results - Results of modulation index calculations (has
    %                       them separated by left and righ sounds)
    %   avg_results - Average results for modulation analysis
    %   sorted_cells - Sorted list of cell indices
    %   all_celltypes - Cell type classifications
    %   selectivity_indexm - Selectivity index matrix
    %   data_type - String specifying the context ('sound' or other)
    %
    % The function iterates over different modulation parameter sets and
    % performs selectivity analysis, saving results to specified directories.
    
    % Define the parameter sets
    param_sets = { 
        struct('mod_threshold', 0.1, 'threshold_single_side', 1, 'savepath', [base 'positive_modulated']),
        struct('mod_threshold', -0.1, 'threshold_single_side', 1, 'savepath', [base 'negative_modulated']),
        struct('mod_threshold', 0.1, 'threshold_single_side', 0, 'savepath', [base 'all_modulated'])
    };
    
    for i = 1:length(param_sets)
        mod_params = param_sets{i};
        if contains(data_type, 'sound')
            mod_params.chosen_mice = 1:25;
        else
            mod_params.chosen_mice = 1:24;
        end
        
        selectivity_params.savepath = mod_params.savepath;
        mkdir(selectivity_params.savepath);
        
        %separate sig cells based on threshold (and single side or not)
        current_sig_cells = get_thresholded_sig_cells(params.info, mod_params, mod_indexm, sig_mod_boot, sorted_cells, all_celltypes, selectivity_params.savepath, 0);
        
        %decide which contexts to use to define significant neurons
        if contains(data_type, 'sound')
            sig_cells = get_significant_neurons(current_sig_cells, mod_indexm, 'union'); %union of active and passive
        else
            sig_cells = get_significant_neurons(current_sig_cells, mod_indexm, 'spont');
        end
        
        %separate neurons into left/right/non selective based on their
        %selecvitiy index (currently using .1 as threshold)
        [selectivity_results_by_dataset, selectivity_results] = analyze_selectivity_pools(selectivity_indexm, mod_indexm, mod_index_results, sig_cells', all_celltypes, params);
        
        % Plot modulation indices for each pool
        generate_selectivity_pool_plots(selectivity_results, selectivity_indexm,mod_index_results, avg_results, all_celltypes, sig_cells, params, selectivity_params.savepath, varargin);
    end
end
