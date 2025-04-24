function flagged_neurons = generate_selectivity_pool_plots(selectivity_results,selectivity_indexm,mod_index_results,avg_results,all_celltypes,sig_cells, params, savepath,varargin)
% scatter of selectivity index across contexts (significant neurons only)
modl_fit = scatter_index_sigcells(sig_cells, all_celltypes, [{selectivity_indexm{:,1}}',{selectivity_indexm{:,2}}'], params.plot_info, savepath, 'Active Selectivity', 'Passive Selectivity');

% histogram of selectivity index (significant neurons only)
[p_val_mod] = histogram_diff_index_sig_cells(sig_cells, all_celltypes,  [{selectivity_indexm{:,1}}',{selectivity_indexm{:,2}}'], params.plot_info, savepath, 'Abs(active) - Abs(passive');

%heatmap separated by selectivity
plot_selectivity_comparison(selectivity_results, savepath); 

%scatter of mod separated by selectivity
plot_selectivity_consistency(selectivity_results, savepath); 

%plots counts of preferred side
flagged_neurons = plot_side_preference(selectivity_results, params, savepath); 

 %scatter plot of modulation index separated by sides and selectivity
scatter_selectivity_vs_modulation(selectivity_results, mod_index_results, savepath);

%average plots separated by selectivity
if nargin > 8
    plot_avg_traces_direction_comparison(avg_results, selectivity_results, savepath,varargin);
    plot_avg_heatmap_direction_comparison(avg_results, selectivity_results, savepath,varargin);
else
    plot_avg_traces_direction_comparison(avg_results, selectivity_results, savepath);
    plot_avg_heatmap_direction_comparison(avg_results, selectivity_results, savepath);
end

