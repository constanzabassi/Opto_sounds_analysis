addpath(genpath('C:\Code\Github\Opto_sounds_analysis'))
% Setup analysis parameters
%includes all datasets being analyzed, frame parameters, mod index
%parameters
params = experiment_config(); 
plot_info = plotting_config(); %plotting params
params.plot_info = plot_info;

%% Pool activity across mice
[all_celltypes, dff_st, deconv_st, stim_info, ...
 mouse_context_tr, deconv_st_interp, alignment_frames] = ...
    pool_activity(params.info.mouse_date, params.info.serverid, params.info.path_string, true, [60,60]);

% Save basic information
% Save variables with consistent paths
filename = fullfile(params.info.savepath, 'data_info');
save(fullfile(filename, 'info.mat'), 'params');
save(fullfile(filename, 'alignment_frames.mat'), 'alignment_frames');
save(fullfile(filename, 'mouse_context_tr.mat'), 'mouse_context_tr');
save(fullfile(filename, 'stim_info.mat'), 'stim_info');

% Process cell types
[num_cells, sorted_cells] = organize_pooled_celltypes(dff_st, all_celltypes); %gives index relative to all datasets
save(fullfile(filename, 'sorted_cells.mat'), 'sorted_cells');

% Separate neural data into contexts (organizes original _st into cell arrays separated by
% context using mouse_context_trials)
%organized context.dff{context,mouse};
[context_data.dff,stim_trials_context,ctrl_trials_context] = separate_structure_2context(dff_st,mouse_context_tr,stim_info);%  context.dff{context,mouse}
[context_data.deconv] = separate_structure_2context(deconv_st,mouse_context_tr,stim_info);%  context.dff{context,mouse}
[context_data.deconv_interp] = separate_structure_2context(deconv_st_interp,mouse_context_tr,stim_info);%  context.dff{context,mouse}

%% Get average responses
% Setup parameters
avg_params = struct(...
    'response_window', 1:122, ...
    'trial_type', 'stim', ...
    'mode', 'separate');

% Get averages
[avg_results_stim ,avg_results_by_dataset_stim,avg_results_ctrl, avg_results_by_dataset_ctrl] = wrapper_trial_averaging(params.info, dff_st,stim_trials_context,ctrl_trials_context, avg_params, [params.info.savepath '/avg/']);
% generate_heatmaps(context_data, sorted_cells, info);

%% Calculate modulation indices
mod_params = params.mod;
mod_params.savepath = fullfile(params.info.savepath, 'mod', mod_params.mod_type, mod_params.mode);

[mod_index_results, sig_mod_boot, mod_indexm] = ...
    wrapper_mod_index_calculation(params.info, dff_st, mod_params.response_range, mod_params.mod_type, mod_params.mode, stim_trials_context, ctrl_trials_context,mod_params.nShuffles, mod_params.savepath);
%% Compare modulation indices across contexts and cell types
mod_params.mod_threshold = .1;% 0 is no threshold applied
mod_params.chosen_mice = [1:24];

%plot % modulated cells per context
sig_mod_boot_thr = plot_pie_thresholded_mod_index(params.info, mod_params, mod_indexm, sig_mod_boot, sorted_cells,mod_params.savepath);
% sig_mod_boot_thr_spont = plot_pie_thresholded_mod_index(info, mod_params, mod_indexm(:,3), sig_mod_boot(:,3), sorted_cells,fullfile(mod_params.savepath,'spont_sig'));

%plot % overlap of modulated cells across contexts!
contexts_to_compare = [1,2]; %[1:3];%[1,2]; %[1,2]; %[1:3];
overlap_labels = {'Active', 'Passive','Both'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %
percent_cells = calculate_sig_overlap_percentages(sig_mod_boot_thr(mod_params.chosen_mice,:), mod_indexm, contexts_to_compare);
plot_sig_overlap_pie(percent_cells, overlap_labels, mod_params.savepath, contexts_to_compare);

% ORGANIZE MODULATION INDICES AND CELL TYPE INDICES ACROSS DATASETS
[context_mod_all, chosen_pyr, chosen_mcherry, chosen_tdtom, celltypes_ids] = ...
    organize_sig_mod_index_contexts_celltypes([1:24], mod_indexm, sig_mod_boot_thr, all_celltypes,plot_info.celltype_names);

%% Make plots of modulation index across contexts/cell types
% Set y-axis limits for the plots.
plot_info.y_lims = [-.4, .4];
params.info.chosen_mice = mod_params.chosen_mice;
% Set labels for plots.
plot_info.plot_labels = {'Stim','Ctrl'}; % Alternative could be {'Left Sounds','Right Sounds'}
plot_info.behavioral_contexts = {'Active','Passive'}; %decide which contexts to plot
overlap_labels = {'Active', 'Passive','Both'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %
params.plot_info = plot_info;

%save directory
save_dir = [mod_params.savepath];% '/spont_sig'];% '/spont_sig']; %[info.savepath '/mod/' mod_params.mod_type '/spont_sig']; % Set directory to save figures.

%generates heatmaps, cdf, box plots, scatter of abs(mod _index)
mod_index_stats = plot_context_comparisons(contexts_to_compare,overlap_labels, mod_indexm, sig_mod_boot, all_celltypes, params,save_dir);

% %scatter plot of modulated cells
% [combined_sig_cells, ~] = union_sig_cells(sig_mod_boot_thr(:,1)', sig_mod_boot_thr(:,2)', mod_indexm);
% modl_fit = scatter_index_sigcells(combined_sig_cells, all_celltypes, [{mod_indexm{:,1}}',{mod_indexm{:,2}}'], plot_info, [], 'Active Mod', 'Passive Mod')

%% Save Results- save your modulation index data.
save(fullfile(save_dir, 'mod_index_data.mat'), 'context_mod_all', 'chosen_pyr', 'chosen_mcherry', 'chosen_tdtom', 'celltypes_ids');
save(fullfile(save_dir, 'mod_index_stats.mat'), 'mod_index_stats');








