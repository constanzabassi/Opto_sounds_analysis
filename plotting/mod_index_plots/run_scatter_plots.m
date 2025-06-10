%% make scatter plots like chengcheng
%0) load data
%opto
opto_sig_mod_boot_thr = load('V:\Connie\results\opto_sound_2025\context\mod\prepost\separate\sig_mod_boot_thr.mat').sig_mod_boot_thr;
opto_sig_mod_boot = load('V:\Connie\results\opto_sound_2025\context\mod\prepost\separate\sig_mod_boot.mat').sig_mod_boot;
opto_mod = load('V:\Connie\results\opto_sound_2025\context\mod\ctrl\separate\mod_indexm.mat').mod_indexm;
opto_mod_prepost = load('V:\Connie\results\opto_sound_2025\context\mod\prepost\separate\mod_indexm.mat').mod_indexm;
opto_mod_num = load('V:\Connie\results\opto_sound_2025\context\mod\ctrl_num\separate\mod_indexm.mat').mod_indexm; %numerator


opto_mod_results = load('V:\Connie\results\opto_sound_2025\context\mod\ctrl\separate\mod_index_results.mat').results;
opto_sig_cells = opto_sig_mod_boot_thr(:,3); %from spontaneous context
opto_average = load('V:\Connie\results\opto_sound_2025\context\avg\trial_averaged_results.mat').avg_results;

%sound
sound_sig_mod_boot_thr = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\sig_mod_boot_thr.mat').sig_mod_boot_thr;
sound_sig_mod_boot = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\sig_mod_boot.mat').sig_mod_boot;
sound_mod = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\mod_indexm.mat').mod_indexm;
sound_mod_num = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound_num\separate\mod_indexm.mat').mod_indexm;

sound_mod_results = load('V:\Connie\results\opto_sound_2025\context\sounds\mod\prepost_sound\separate\mod_index_results.mat').results;
[sound_sig_cells, ~] = union_sig_cells(sound_sig_mod_boot_thr(:,1)', sound_sig_mod_boot_thr(:,2)', sound_mod);
sound_average = load('V:\Connie\results\opto_sound_2025\context\sounds\avg\trial_averaged_results_sounds.mat').avg_results_sounds;
sound_stim_average = load('V:\Connie\results\opto_sound_2025\context\sounds\avg\trial_averaged_results_sounds_stim.mat').avg_results;
%%
%1) use modulation index and plot that across active and passive
save_dir = 'V:\Connie\results\opto_sound_2025\context\scatter_plots'
modl_fit = scatter_index_sigcells(opto_sig_cells, all_celltypes, [{opto_mod{:,1}}',{opto_mod{:,2}}'], plot_info, save_dir, 'Active Stim Mod', 'Passive Stim Mod');
modl_fit = scatter_index_sigcells(opto_sig_cells, all_celltypes, [{opto_mod_num{:,1}}',{opto_mod_num{:,2}}'], plot_info, save_dir, 'Active (stim+sound) - sound', 'Passive (stim+sound) - sound');

modl_fit = scatter_index_sigcells(sound_sig_cells, all_celltypes, [{sound_mod{:,1}}',{sound_mod{:,2}}'], plot_info, save_dir, 'Active Sound Mod', 'Passive Sound Mod');
modl_fit = scatter_index_sigcells(sound_sig_cells, all_celltypes, [{sound_mod_num{:,1}}',{sound_mod_num{:,2}}'], plot_info, save_dir, 'Active post-pre (sound)', 'Passive post-pre (sound)');

