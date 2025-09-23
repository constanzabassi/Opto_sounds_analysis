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
[num_cells, sorted_cells] = organize_pooled_celltypes(dff_st, all_celltypes); %gives index relative to all datasets
[context_data.dff,stim_trials_context,ctrl_trials_context] = separate_structure_2context(dff_st,mouse_context_tr,stim_info);%  context.dff{context,mouse}
[context_data.deconv] = separate_structure_2context(deconv_st,mouse_context_tr,stim_info);%  context.dff{context,mouse}
[context_data.deconv_interp] = separate_structure_2context(deconv_st_interp,mouse_context_tr,stim_info);%  context.dff{context,mouse}
%% Calculate engagement index in the pre stimulus period
mod_params = params.mod;
mod_params.mod_type = 'pre_engagement';
mod_params.mode = 'simple'; %
mod_params.savepath = fullfile(params.info.savepath, 'mod', mod_params.mod_type, mod_params.mode)

[mod_index_results, sig_mod_boot, mod_indexm] = ...
    wrapper_engagement_index_calculation(params.info, dff_st, mod_params.response_range, mod_params.mod_type, mod_params.mode, stim_trials_context, ctrl_trials_context,mod_params.nShuffles, mod_params.savepath);
%% get thresholded significant cells
mod_params.mod_threshold = .1;% 0 is no threshold applied
mod_params.chosen_mice = [1:25];

%using previously calculated mod index from prepost (looking at spont to
%find the thresholded cells!)
load('W:\Connie\results\Bassi2025\fig3\mod\pre_engagement\simple\sig_mod_boot.mat');
load('W:\Connie\results\Bassi2025\fig3\mod\pre_engagement\simple\mod_indexm.mat');
%plot % modulated cells per context
mod_params.threshold_single_side = 0;
sig_mod_boot_thr = get_thresholded_sig_cells(params.info, mod_params, mod_indexm', sig_mod_boot', sorted_cells, all_celltypes, mod_params.savepath,0);
%%
contexts_to_compare = [1]; %[1:3];%[1,2]; %[1,2]; %[1:3];
overlap_labels = {'Act - Pass'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %

% ORGANIZE MODULATION INDICES AND CELL TYPE INDICES ACROSS DATASETS
[context_mod_all, chosen_pyr, chosen_mcherry, chosen_tdtom, celltypes_ids] = ...
    organize_sig_mod_index_contexts_celltypes([1:25], mod_indexm', sig_mod_boot_thr, all_celltypes,plot_info.celltype_names);

[~, ~, ~, ~, celltypes_ids] = ...
    organize_sig_mod_index_contexts_celltypes([1:24], mod_indexm', sig_mod_boot_thr, all_celltypes,plot_info.celltype_names);

%% make plots!
plot_info.y_lims = [-.3, .3];
params.plot_info = plot_info;
params.info.chosen_mice = [1:24];
% Set labels for plots.
plot_info.plot_labels = {'Stim','Ctrl'}; % Alternative could be {'Left Sounds','Right Sounds'}
plot_info.behavioral_contexts = {'Engagement Index'}; %decide which contexts to plot
overlap_labels = {'Active', 'Passive','Both'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %
params.plot_info = plot_info;

%generate average plots
savepath = 'W:\Connie\results\Bassi2025\fig3\pre_engagement\celltype_traces\';
wrapper_avg_cell_type_traces_engagement(context_data.dff,all_celltypes,mod_indexm,sig_mod_boot,mod_params,savepath,'engagement_dff',plot_info, plot_info.celltype_names,plot_info.colors_celltypes_3contexts);

%generates heatmaps, cdf, box plots, scatter of abs(mod _index)
mod_index_stats_datasets = generate_engagement_index_plots_datasets(params.info.chosen_mice, mod_indexm',  sig_mod_boot_thr, all_celltypes, params, mod_params.savepath, celltypes_ids,plot_info.y_lims);
save(fullfile(save_dir, 'mod_index_stats_datasets.mat'), 'mod_index_stats_datasets');

%% compare engagement indices in opto or sound neurons
[sound,opto,sorted_cells,all_celltypes,context_data,ctrl_trials_context,stim_trials_context] = load_processed_opto_sound_data(params,{'separate','separate'});

% set up colors and different pools of cells
%within context use below
% [pooled_cell_types,plot_info.celltype_names,plot_info.colors_celltypes] = organize_functional_groups(all_celltypes, sound.sig_mod_boot_thr, opto.sig_mod_boot_thr_ctrl, opto.mod(1:24,:), {'unmodulated','both','opto','sound'},[1:24],plot_info, 2);

% % significance based on union active/passive for sound, spontaneous for stim
[pooled_cell_types,plot_info.functional_names,plot_info.functional_colors] = organize_functional_groups(all_celltypes, sound.sig_cells, opto.sig_cells, opto.mod(1:24,:), {'sound','opto','both','unmodulated'},[1:24],plot_info, 1);

savepath = 'W:\Connie\results\Bassi2025\fig3\pre_engagement\functional_celltype_traces\';
mod_params.chosen_mice = [1:24]; %1 less for opto control
wrapper_avg_cell_type_traces_engagement(context_data.dff,pooled_cell_types,mod_indexm,sig_mod_boot,mod_params,savepath,'engagement_dff',plot_info,plot_info.functional_names,plot_info.colors_pooled_3contexts); %repelem(plot_info.functional_colors, 3, 1)

