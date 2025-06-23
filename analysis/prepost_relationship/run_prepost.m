
addpath(genpath('C:\Code\Github\Opto_sounds_analysis'))
% Setup analysis parameters
%includes all datasets being analyzed, frame parameters, mod index
%parameters
params = experiment_config(); 
plot_info = plotting_config(); %plotting params
params.plot_info = plot_info;

%% load mod indices/ significant neurons/ selectivity indices

%opto
opto_sig_mod_boot_thr = load('V:\Connie\results\opto_sound_2025\context\mod\prepost\separate\sig_mod_boot_thr.mat').sig_mod_boot_thr;
opto_sig_mod_boot = load('V:\Connie\results\opto_sound_2025\context\mod\prepost\separate\sig_mod_boot.mat').sig_mod_boot;
opto_mod = load('V:\Connie\results\opto_sound_2025\context\mod\ctrl\separate\mod_indexm.mat').mod_indexm;
opto_mod_prepost = load('V:\Connie\results\opto_sound_2025\context\mod\prepost\separate\mod_indexm.mat').mod_indexm;

opto_mod_results = load('V:\Connie\results\opto_sound_2025\context\mod\ctrl\separate\mod_index_results.mat').results;
opto_sig_cells = opto_sig_mod_boot_thr(:,3); %from spontaneous context
opto_average = load('V:\Connie\results\opto_sound_2025\context\avg\trial_averaged_results.mat').avg_results;

%sound
sound_sig_mod_boot_thr = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\sig_mod_boot_thr.mat').sig_mod_boot_thr;
sound_sig_mod_boot = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\sig_mod_boot.mat').sig_mod_boot;
sound_mod = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\mod_indexm.mat').mod_indexm;
sound_mod_results = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\mod_index_results.mat').results;
[sound_sig_cells, ~] = union_sig_cells(sound_sig_mod_boot_thr(:,1)', sound_sig_mod_boot_thr(:,2)', sound_mod);
sound_average = load('V:\Connie\results\opto_sound_2025\context\sounds\avg\trial_averaged_results_sounds.mat').avg_results_sounds;
sound_stim_average = load('V:\Connie\results\opto_sound_2025\context\sounds\avg\trial_averaged_results_sounds_stim.mat').avg_results;

%previously calculated selecvitiy
selectivity_mode = 'prepost_ctrl';
selectivity_base = ['V:\Connie\results\opto_sound_2025\context\sounds\selectivity\' selectivity_mode];
selectivity_indexm = load([selectivity_base '\mod_indexm.mat']).mod_indexm;
selectivity_results = load([selectivity_base '\mod_index_results.mat']).results; %has significant neurons too!

% selectivity_indexm = load('V:\Connie\results\opto_sound_2025\context\sounds\selectivity\prepost_ctrl\mod_indexm.mat').mod_indexm;

load('V:\Connie\results\opto_sound_2025\context\data_info\sorted_cells.mat');
load('V:\Connie\results\opto_sound_2025\context\data_info\all_celltypes.mat');
%load data and trials
load('V:\Connie\results\opto_sound_2025\context\data_info\context_data.mat');
load('V:\Connie\results\opto_sound_2025\context\data_info\ctrl_trials_context.mat');
load('V:\Connie\results\opto_sound_2025\context\data_info\stim_trials_context.mat');

%%

frame_window = 50:59;%63:92;

plot_info.colors_celltypes = [0.3,0.2,0.6 ; 1,0.7,0; 0.3,0.8,1]; %0.9/0.6/0.2 or 1,0.7,0
opto_only_sig_cells = setdiff_sig_cells(opto_sig_cells', sound_sig_cells,opto_mod);
sound_only_sig_cells = setdiff_sig_cells(sound_sig_cells(1,1:24),opto_sig_cells', opto_mod);
opto_sound_sig_cells = intersect_sig_cells(opto_sig_cells', sound_sig_cells,opto_mod);
all_cells =[cellfun(@(x) x.pyr_cells,all_celltypes,'UniformOutput',false);cellfun(@(x) x.som_cells,all_celltypes,'UniformOutput',false);cellfun(@(x) x.pv_cells,all_celltypes,'UniformOutput',false)];
num_cells = cellfun(@length, all_cells );

pooled_cell_types = {};
for dataset = 1:24
    current_sig_cells = [sound_only_sig_cells{dataset}, opto_only_sig_cells{dataset}, opto_sound_sig_cells{dataset}];
%     pooled_cell_types{dataset}.unmodulated = setdiff(1:sum(num_cells(:,dataset)),current_sig_cells);
    pooled_cell_types{dataset}.sound = sound_only_sig_cells{dataset};
    pooled_cell_types{dataset}.opto = opto_only_sig_cells{dataset};
    pooled_cell_types{dataset}.both = opto_sound_sig_cells{dataset};

end

%% difference of post and pre
plot_info.celltype_names = {'Sound','Opto','Both'};
plot_info.y_lims = [-.2, .4];
% Set labels for plots.
plot_info.behavioral_contexts = {'Active','Passive'}; %decide which contexts to plot
params.plot_info = plot_info;

[dff_response,~] = unpack_context_mouse_celltypes(context_data.dff,[],all_celltypes,[1:25]); %context_data.deconv_interp
[dff_trial_cel_mouse_pre,dff_context_celltype] = calc_avg_rate_across_context_celltype_choosetrials(dff_response,50:59,stim_trials_context, ctrl_trials_context);
[dff_trial_cel_mouse_post,dff_context_celltype2] = calc_avg_rate_across_context_celltype_choosetrials(dff_response,63:92,stim_trials_context, ctrl_trials_context);

% if wanted to do deconvolved
% [deconv_response,~] = unpack_context_mouse_celltypes(context_data.deconv_interp,[],all_celltypes,[1:25]); %context_data.deconv_interp
% [dff_trial_cel_mouse_pre,dff_context_celltype] = calc_spike_rate_across_context_celltype_choosetrials(deconv_response,50:59,stim_trials_context, ctrl_trials_context);
% [dff_trial_cel_mouse_post,dff_context_celltype2] = calc_spike_rate_across_context_celltype_choosetrials(deconv_response,63:92,stim_trials_context, ctrl_trials_context);



[dff_context_matrix_pre,dff_context_matrix_ctrl_pre] = into_mod_structure(dff_trial_cel_mouse_pre,all_celltypes); %same as mod structure for easy plotting!
[dff_context_matrix_post,dff_context_matrix_ctrl_post] = into_mod_structure(dff_trial_cel_mouse_post,all_celltypes); %same as mod structure for easy plotting!

% define pre and post periods
[diff_stim, diff_pre_stim, diff_pre_ctrl] = calculate_avg_differences(dff_context_matrix_pre,dff_context_matrix_ctrl_pre,dff_context_matrix_post,dff_context_matrix_ctrl_post);

%make scatter plots and save them!
current_save_dir = ''; %'V:\Connie\results\opto_sound_2025\context\mod_index_specified_cells\differences_pre_post\dff';
modl_fit = scatter_index_sigcells_histogram([], pooled_cell_types, [{diff_pre_stim{:,1}}',{diff_stim{:,1}}'], plot_info, current_save_dir, 'Pre Diff (active - passive)', 'Active Post (stim+sound - sound)',[-.6,2]);
modl_fit = scatter_index_sigcells_histogram([], pooled_cell_types, [{diff_pre_stim{:,1}}',{diff_stim{:,2}}'], plot_info, current_save_dir, 'Pre Diff (active - passive)', 'Passive Post (stim+sound - sound)',[-.6,2]);

modl_fit = scatter_index_sigcells_histogram([], pooled_cell_types, [{diff_pre_ctrl{:,1}}',{dff_context_matrix_ctrl_post{:,1}}'], plot_info, current_save_dir, 'Pre Diff (active - passive)', 'Active Post (sound)',[-.6,2]);
modl_fit = scatter_index_sigcells_histogram([], pooled_cell_types, [{diff_pre_ctrl{:,1}}',{dff_context_matrix_ctrl_post{:,2}}'], plot_info, current_save_dir, 'Pre Diff (active - passive)', 'Passive Post (sound)',[-.6,2]);

[preavg_index_by_dataset,~] = unpack_modindexm(dff_context_matrix_pre,[],pooled_cell_types,[1:24]);
preavg_stats_celltypes_dataset = plot_connected_abs_mod_by_mouse(current_save_dir, preavg_index_by_dataset, [1:24],...
          params.plot_info, [0,.3],0,'Pre Mean (ΔF/F)');

%across celltypes
plot_info.colors_celltypes = [0.37 0.75 0.49 %light green
                            0.17 0.35 0.8  %blue
                            0.82 0.04 0.04]; % red  
[preavg_index_by_dataset_ct,~] = unpack_modindexm(dff_context_matrix_pre,[],all_celltypes,[1:24]);
preavg_stats_celltypes_dataset_ct = plot_connected_abs_mod_by_mouse([current_save_dir '/celltypes'], preavg_index_by_dataset_ct, [1:24],...
          plot_info, [0,.3],0,'Pre Mean (ΔF/F)');