addpath(genpath('C:\Code\Github\Opto_sounds_analysis'))
% Setup analysis parameters
%includes all datasets being analyzed, frame parameters, mod index
%parameters
params = experiment_config(); 
plot_info = plotting_config(); %plotting params
params.plot_info = plot_info;
%% load mod indices/ significant neurons
[sound,opto,sorted_cells,all_celltypes,context_data,ctrl_trials_context,stim_trials_context] = load_processed_opto_sound_data(params,{'separate','separate'});
selectivity_mode = 'prepost_ctrl'; %or prepost_ctrl or ctrl (just difference between left and right sounds)
selectivity_base = ['V:\Connie\results\opto_sound_2025\context\sounds\selectivity\' selectivity_mode];
selectivity_indexm = load([selectivity_base '\mod_indexm.mat']).mod_indexm;
selectivity_results = load([selectivity_base '\mod_index_results.mat']).results; %has significant neurons too!

load('V:\Connie\results\opto_sound_2025\context\data_info\sorted_cells.mat');
%% Analyze modulation indices by selectivity pools
%Sounds first
mod_indexm = sound.mod;
sig_mod_boot = sound.sig_mod_boot;
mod_index_results = sound.results;
avg_results = sound.avg;
data_type = 'sounds';
params.selectivity_sounds.selectivity_sig_mode = 'union'; %'union' or 'intersect'
base = ['W:\Connie\results\Bassi2025\fig3\selectivity_pools\' data_type '\' selectivity_mode '\' params.selectivity_sounds.selectivity_sig_mode '\'];% 'V:\Connie\results\opto_sound_2025\context\sounds\selectivity\negative';
mkdir(base);

[selectivity_pool_results_by_dataset, selectivity_pool_results] = wrapper_selecitivity_pool_analysis(base, params, mod_indexm,[], sig_mod_boot, mod_index_results,selectivity_results, avg_results, sorted_cells, all_celltypes, selectivity_indexm, data_type,[.1,.4],'Avg. Sound (ΔF/F)');
% results organized as selectivity{i} where i is +, -, or all modulatedneurons
save('selectivity_pool_results','selectivity_pool_results_by_dataset', 'selectivity_pool_results');
save('selectivity_pool_results','selectivity_pool_results');
wrapper_selecitivity_pool_analysis(base, params, mod_indexm,[], sig_mod_boot, mod_index_results,selectivity_results, sound.stim_avg, sorted_cells, all_celltypes, selectivity_indexm, data_type,[.1,.4],'Avg. Sound+Stim (ΔF/F)');
%% Analyze modulation indices by selectivity pools
%OPTO!
mod_indexm = opto.mod;
sig_mod_boot = opto.sig_mod_boot;
mod_index_results = opto.results;
avg_results = find_average_difference(opto.avg, sound.avg); %opto_average;
data_type = 'opto';
params.selectivity_sounds.selectivity_sig_mode = 'union'; %'union' or 'intersect'
base = ['W:\Connie\results\Bassi2025\fig3\selectivity_pools\' data_type '\' selectivity_mode '\' params.selectivity_sounds.selectivity_sig_mode '\'];% 'V:\Connie\results\opto_sound_2025\context\sounds\selectivity\negative';
mkdir(base);

[selectivity_pool_results_by_dataset, selectivity_pool_results] = wrapper_selecitivity_pool_analysis(base, params, mod_indexm,opto.mod_prepost, sig_mod_boot, mod_index_results, selectivity_results, avg_results, sorted_cells, all_celltypes, selectivity_indexm, data_type,[-.2,.5],'Difference in ΔF/F');
save('selectivity_pool_results','selectivity_pool_results_by_dataset', 'selectivity_pool_results');

avg_results = opto.avg;
wrapper_selecitivity_pool_analysis(base, params, mod_indexm, opto.mod_prepost, sig_mod_boot, mod_index_results, selectivity_results, avg_results, sorted_cells, all_celltypes, selectivity_indexm, data_type,[-.2,.5],'Avg. Stim+Sound ΔF/F');

avg_results = sound.avg;
wrapper_selecitivity_pool_analysis(base, params, mod_indexm, opto.mod_prepost, sig_mod_boot, mod_index_results, selectivity_results, avg_results, sorted_cells, all_celltypes, selectivity_indexm, data_type,[-.2,.5],'Avg. Sound ΔF/F');

%% compare the overlap of opto neurons and sound neurons

contexts_to_compare = [1,2]; %[1:3];%[1,2]; %[1,2]; %[1:3];
overlap_labels = {'Photostim Only', 'Sound Only','Photostim & Sound','Not Modulated'}; 
save_dir = ['W:\Connie\results\Bassi2025\fig3\'];
[percent_cells, percent_cells_per_dataset] = calculate_sig1_vs_sig2_overlap(opto.sig_cells(1:24,:), sound.sig_mod_boot_thr, opto.mod, contexts_to_compare);
plot_sig_overlap_pie(mean(percent_cells_per_dataset), overlap_labels, save_dir, contexts_to_compare);

%% compare sound vs photostim selectivity indices
[mod_index_by_dataset,~] = unpack_modindexm(selectivity_indexm,sound.sig_cells,all_celltypes,[1:24]);
[mod_index_by_dataset_opto,~] = unpack_modindexm(selectivity_indexm,opto.sig_cells,all_celltypes,[1:24]);

 %by dataset
 save_dir = ['W:\Connie\results\Bassi2025\fig3\selectivity_pools\'];
% all_stats.abs_mod_stats_celltypes_dataset = plot_connected_abs_mod_by_mouse(save_dir, mod_index_by_dataset, [1:24],...
%   params.plot_info, [0,.5]);
abs_logic = 1;
avg_contexts = 1;
mod_stats = compare_mod_indices_across_conditions(save_dir, ...
mod_index_by_dataset, mod_index_by_dataset_opto, ...
[1:24], [1:24], params.plot_info, [0,.3],abs_logic,avg_contexts);
