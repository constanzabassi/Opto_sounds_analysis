contexts_to_compare = [1,2];

params = experiment_config(); 
plot_info = plotting_config(); %plotting params
params.plot_info = plot_info;
load('V:\Connie\results\opto_sound_2025\context\data_info\all_celltypes.mat')
%% Sound Index Plots
mod_params = params.mod_sounds;
mod_params.mod_threshold = .1;% 0 is no threshold applied
mod_params.chosen_mice = [1:25];
% 1) load data
load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\mod_indexm.mat');
load('V:\Connie\results\opto_sound_2025\context\sounds\mod\separate\sig_mod_boot_thr.mat')
%%%% sig cells %%%%%%%%%%%
plot_info.y_lims = [-.4, .4];
% Set labels for plots.
plot_info.plot_labels = {'Sounds','Sounds'}; % Alternative could be {'Left Sounds','Right Sounds'}
plot_info.behavioral_contexts = {'Active','Passive'}; %decide which contexts to plot
overlap_labels = {'Active', 'Passive','Both'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %
params.plot_info = plot_info;
params.info.chosen_mice = [1:25];
params.string = 'Sounds';

%%%% sig cells %%%%%%%%%%%
savepath = 'W:\Connie\results\Bassi2025\fig3\sounds\mod\prepost_sound\separate\sig_neurons';
[combined_sig_cells, ~] = union_sig_cells(sig_mod_boot_thr(:,1)', sig_mod_boot_thr(:,2)', mod_indexm);
%single neurons
plot_info.y_lims = [-.4, .4];params.plot_info = plot_info;
mod_index_stats = plot_context_comparisons(contexts_to_compare,overlap_labels, mod_indexm, combined_sig_cells, all_celltypes, params, save_dir);
%datasets
plot_info.y_lims = [-.2, .3];params.plot_info = plot_info;
mod_index_stats_datasets = generate_mod_index_plots_datasets(params.info.chosen_mice, mod_indexm, combined_sig_cells, all_celltypes, params, save_dir);
%%% STATS TABLES
% table_fig3 = make_stats_tables_mod_index(mod_index_stats, mod_index_stats_datasets, save_dir);

%%%% all cells %%%%%%%%%%%
%save directory
save_dir = 'W:\Connie\results\Bassi2025\fig3\sounds\mod\prepost_sound\separate\';% '/spont_sig'];% '/spont_sig']; %[info.savepath '/mod/' mod_params.mod_type '/spont_sig']; % Set directory to save figures.

%single neurons
plot_info.y_lims = [-.4, .4];params.plot_info = plot_info;
mod_index_stats = plot_context_comparisons(contexts_to_compare,overlap_labels, mod_indexm, [], all_celltypes, params, save_dir);
plot_info.y_lims = [-.2, .20];params.plot_info = plot_info;
%datasets
mod_index_stats_datasets = generate_mod_index_plots_datasets(params.info.chosen_mice, mod_indexm, [], all_celltypes, params, save_dir);
%%% STATS TABLES
% table_fig3 = make_stats_tables_mod_index(mod_index_stats, mod_index_stats_datasets, save_dir);
%% Photostim Index plots
mod_params = params.mod;
mod_params.mod_threshold = .1;% 0 is no threshold applied
mod_params.chosen_mice = [1:24];
params.string = 'opto';
% 1) load data
load('V:\Connie\results\opto_sound_2025\context\mod\prepost\separate\sig_mod_boot_thr.mat')% sig neurons based on pre post spont
load('V:\Connie\results\opto_sound_2025\context\mod\ctrl\separate\mod_indexm.mat')
% Set y-axis limits for the plots.
plot_info.y_lims = [-.4, .4];
params.info.chosen_mice = mod_params.chosen_mice;
% Set labels for plots.
plot_info.plot_labels = {'Stim','Ctrl'}; % Alternative could be {'Left Sounds','Right Sounds'}
plot_info.behavioral_contexts = {'Active','Passive'}; %decide which contexts to plot
overlap_labels = {'Active', 'Passive','Both'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %
params.plot_info = plot_info;

%%%% sig cells %%%%%%%%%%%
savepath = 'W:\Connie\results\Bassi2025\fig3\mod\ctrl\separate\sig_neurons';
%single neurons
mod_index_stats = plot_context_comparisons(contexts_to_compare,overlap_labels, mod_indexm, sig_mod_boot_thr(:,3), all_celltypes, params,savepath); %single neurons
%datasets
plot_info.y_lims = [-.2, .4]; params.plot_info = plot_info;
mod_index_stats_datasets = generate_mod_index_plots_datasets(params.info.chosen_mice, mod_indexm,  sig_mod_boot_thr(:,3)', all_celltypes, params,savepath);
%%% STATS TABLES
% table_fig3 = make_stats_tables_mod_index(mod_index_stats, mod_index_stats_datasets, save_dir);

%%%% all cells %%%%%%%%%%%
plot_info.y_lims = [-.4, .4];params.plot_info = plot_info;
%save directory
savepath = 'W:\Connie\results\Bassi2025\fig3\mod\ctrl\separate';%

%single neurons
mod_index_stats = plot_context_comparisons(contexts_to_compare,overlap_labels, mod_indexm, [], all_celltypes, params,savepath);%single neurons
plot_info.y_lims = [-.2, .2];params.plot_info = plot_info;
%datasets
mod_index_stats_datasets = generate_mod_index_plots_datasets(params.info.chosen_mice, mod_indexm,  [], all_celltypes, params,savepath);

%%% STATS TABLES
% table_fig3 = make_stats_tables_mod_index(mod_index_stats, mod_index_stats_datasets, save_dir);

%% Engagement Index Plots
% 1) load data

% sig cells


% all cells
