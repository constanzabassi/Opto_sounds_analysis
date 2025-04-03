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

%previously calculated selecvitiy
% selectivity_indexm = load('V:\Connie\results\opto_sound_2025\context\sounds\selectivity\mod_indexm.mat').mod_indexm;
selectivity_indexm = load('V:\Connie\results\opto_sound_2025\context\sounds\selectivity\prepost_ctrl\mod_indexm.mat').mod_indexm;

load('V:\Connie\results\opto_sound_2025\context\data_info\sorted_cells.mat');
load('V:\Connie\results\opto_sound_2025\context\data_info\all_celltypes.mat');
%% Analyze modulation indices by selectivity pools
%Sounds first
mod_indexm = sound_mod;
sig_mod_boot = sound_sig_mod_boot;
mod_index_results = sound_mod_results;
avg_results = sound_average;
data_type = 'sounds';
base = ['V:\Connie\results\opto_sound_2025\context\selectivity_pools\' data_type '\prepost_ctrl\'];% 'V:\Connie\results\opto_sound_2025\context\sounds\selectivity\negative';
mkdir(base);

wrapper_selecitivity_pool_analysis(base, params, mod_indexm, sig_mod_boot, mod_index_results, avg_results, sorted_cells, all_celltypes, selectivity_indexm, data_type,[.1,.3],'Response (ΔF/F)');

%% opto second
mod_indexm = opto_mod;
sig_mod_boot = opto_sig_mod_boot;
mod_index_results = opto_mod_results;
avg_results = find_average_difference(opto_average, sound_average); %opto_average;
data_type = 'opto';
base = ['V:\Connie\results\opto_sound_2025\context\selectivity_pools\' data_type '\prepost_ctrl\'];% 'V:\Connie\results\opto_sound_2025\context\sounds\selectivity\negative';
mkdir(base);

wrapper_selecitivity_pool_analysis(base, params, mod_indexm, sig_mod_boot, mod_index_results, avg_results, sorted_cells, all_celltypes, selectivity_indexm, data_type,[-.1,.2],'Difference in ΔF/F');

%%
% mod_params.mod_threshold = .1;% 0 is no threshold applied
% mod_params.threshold_single_side = 1;% 0 is no threshold applied
% mod_params.chosen_mice = 1:25;
% selectivity_params.savepath = 'V:\Connie\results\opto_sound_2025\context\sounds\selectivity\prepost_ctrl\sound_opto_cells';
% mkdir(selectivity_params.savepath)
