addpath(genpath('C:\Code\Github\Opto_sounds_analysis'))
% Setup analysis parameters
%includes all datasets being analyzed, frame parameters, mod index
%parameters
params = experiment_config(); 
plot_info = plotting_config(); %plotting params
params.plot_info = plot_info;

%% Compile sound and opto trials onsets and find virmen information from imaging structure for those trials
%get sound info (a bit slow so I saved inside:
%'V:/Connie/results/opto_sound_2025/context' loaded in pooled_activity_sounds

[sound_onsets_all, alignment_frames_all, control_output_all, opto_output_all,sound_only_all, loc_trial, all_trial_info_sounds] = compile_trial_data_stim_sound(params.info,{'active','passive'},params.info.savepath);
%% Pool activity across mice
% neural data will be sound + ctrl or sound only trials!!!
[all_celltypes, sound_data, sound_context_data]  = ...
    pool_activity_sounds(params.info.mouse_date, params.info.serverid, [60,60]);

% Save basic information
% Save variables with consistent paths
% filename = fullfile(params.info.savepath, 'data_info');
% save(fullfile(filename, 'info.mat'), 'params');
% save(fullfile(filename, 'alignment_frames.mat'), 'alignment_frames');
% save(fullfile(filename, 'mouse_context_tr.mat'), 'mouse_context_tr');
% save(fullfile(filename, 'stim_info.mat'), 'stim_info');

% Process cell types
[num_cells, sorted_cells] = organize_pooled_celltypes(sound_data.active.dff_st, all_celltypes); %gives index relative to all datasets
% save(fullfile(filename, 'sorted_cells.mat'), 'sorted_cells');

% Separate neural data into contexts (organizes original _st into cell arrays separated by
% context using mouse_context_trials)
%organized context_data.dff{context,mouse};
[context_data.dff,stim_trials_context,ctrl_trials_context] = organize_2context(sound_data.active.dff_st,sound_data.passive.dff_st);
[stim_info_combined,dff_st_combined] = combine_stim_info_dff_st(sound_context_data.active, sound_context_data.passive, sound_data.active.dff_st,sound_data.passive.dff_st);

[context_data.deconv,~,~] = organize_2context(sound_data.active.deconv_st_interp,sound_data.passive.deconv_st_interp);
[~,deconv_st_combined] = combine_stim_info_dff_st(sound_context_data.active, sound_context_data.passive, sound_data.active.deconv_st_interp,sound_data.passive.deconv_st_interp);

%% Get average responses
% Setup parameters
avg_params = struct(...
    'response_window', 1:122, ...
    'trial_type', 'sounds', ...
    'mode', 'separate');

% Get averages
[avg_results_sounds_stim ,avg_results_by_dataset_sounds_stim,avg_results_sounds, avg_results_by_dataset_sounds] = wrapper_trial_averaging(params.info, dff_st_combined,stim_trials_context,ctrl_trials_context, avg_params, []);
% generate_heatmaps(context_data, sorted_cells, info);

%% Calculate modulation indices
mod_params = params.mod_sounds; %use 'prespose'/'separate'?
mod_params.savepath = fullfile(params.info.savepath_sounds, 'mod', mod_params.mod_type, mod_params.mode);

[mod_index_results, sig_mod_boot, mod_indexm] = ...
    wrapper_mod_index_calculation(params.info, dff_st_combined, mod_params.response_range, mod_params.mod_type, mod_params.mode, stim_trials_context, ctrl_trials_context,mod_params.nShuffles, mod_params.simple_or_not, mod_params.savepath);
%% Compare modulation indices across contexts and cell types
mod_params.mod_threshold = .1;% 0 is no threshold applied
mod_params.chosen_mice = [1:25];

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
    organize_sig_mod_index_contexts_celltypes([1:25], mod_indexm, sig_mod_boot_thr, all_celltypes,plot_info.celltype_names);

%% Make plots of modulation index across contexts/cell types
% Set y-axis limits for the plots.
plot_info.y_lims = [-.4, .4];
% Set labels for plots.
plot_info.plot_labels = {'Left Sounds','Right Sounds'}; % Alternative could be {'Left Sounds','Right Sounds'}
plot_info.behavioral_contexts = {'Active','Passive'}; %decide which contexts to plot
overlap_labels = {'Active', 'Passive','Both'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %
params.plot_info = plot_info;
params.info.chosen_mice = [1:25];

%save directory
save_dir = [mod_params.savepath];% '/spont_sig'];% '/spont_sig']; %[info.savepath '/mod/' mod_params.mod_type '/spont_sig']; % Set directory to save figures.

%generates heatmaps, cdf, box plots, scatter of abs(mod _index)
mod_index_stats = plot_context_comparisons(contexts_to_compare,overlap_labels, mod_indexm, sig_mod_boot, all_celltypes, params, save_dir);


%% Save Results- save your modulation index data.
save(fullfile(save_dir, 'mod_index_data.mat'), 'context_mod_all', 'chosen_pyr', 'chosen_mcherry', 'chosen_tdtom', 'celltypes_ids');

%% Now compute selectivity Indices?
selectivity_params = params.selectivity_sounds; 
selectivity_params.savepath = fullfile(params.info.savepath_sounds, 'selectivity');

[selectivity_index_results, selectivity_sig_mod_boot, selectivity_indexm] = ...
    wrapper_mod_index_calculation(params.info, dff_st_combined, selectivity_params.response_range, selectivity_params.mod_type, selectivity_params.mode, stim_trials_context, ctrl_trials_context,selectivity_params.nShuffles, selectivity_params.simple_or_not, selectivity_params.savepath);


[combined_sig_cells, ~] = union_sig_cells(sig_mod_boot_thr(:,1)', sig_mod_boot_thr(:,2)', mod_indexm);
modl_fit = scatter_index_sigcells(combined_sig_cells, all_celltypes, [{selectivity_indexm{:,1}}',{selectivity_indexm{:,2}}'], plot_info, selectivity_params.savepath, 'Active Selectivity', 'Passive Selectivity')
[p_val_mod] = histogram_diff_index_sig_cells(combined_sig_cells, all_celltypes,  [{selectivity_indexm{:,1}}',{selectivity_indexm{:,2}}'], plot_info, selectivity_params.savepath, 'Abs(active) - Abs(passive')

%% Analyze modulation indices by selectivity pools
mod_params.mod_threshold = .1;% 0 is no threshold applied
sig_mod_boot_thr = plot_pie_thresholded_mod_index(params.info, mod_params, mod_indexm, sig_mod_boot, sorted_cells,mod_params.savepath);

[selectivity_results_by_dataset,selectivity_results] = analyze_selectivity_pools(selectivity_indexm, ...
    mod_indexm,mod_index_results, sig_mod_boot_thr, all_celltypes, params);

% Save selectivity results
save(fullfile(selectivity_params.savepath, 'selectivity_results.mat'), 'selectivity_results');

% Plot modulation indices for each pool
plot_selectivity_comparison(selectivity_results, selectivity_params.savepath); %heatmap
plot_selectivity_consistency(selectivity_results, selectivity_params.savepath); %scatter of mod separated by selectivity
plot_side_preference(selectivity_results, params); %plots counts of preferred side
scatter_selectivity_vs_modulation(selectivity_indexm, mod_index_results); %scatter plot of modulation index separated by sides and selectivity