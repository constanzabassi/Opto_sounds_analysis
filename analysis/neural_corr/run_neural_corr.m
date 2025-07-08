corr_params.corr_bins = 3; %how many frames to bin
corr_params.nonoverlap_bins = 1;
corr_params.type = 'stim'; %stim or sound
corr_params.field_name_neural = 'stim'; %stim or ctrl
corr_params.calc_type = 'snr'; %corr or snr
save_data_directory =['V:\Connie\results\opto_sound_2025\context\snr\opto'];


[mouse_corr_stats, trial_indices_right_all, trial_indices_left_all] = neural_paircorr_acrossmice([1:24], params.info.mouse_date, params.info.serverid, dff_st,combined_sig_cells, [], 63:92,stim_trials_context, ctrl_trials_context,  corr_params,save_data_directory);

%make cdf plot
p_values = cdf_corr_across_contexts(mouse_corr_stats, params.plot_info, save_data_directory,[-1:0.05:1]);
%%
corr_params.corr_bins = 3; %how many frames to bin
corr_params.nonoverlap_bins = 1;
corr_params.type = 'stim'; %stim or sound
corr_params.field_name_neural = 'stim'; %stim or ctrl
corr_params.calc_type = 'noise'; %corr or snr
save_data_directory =['V:\Connie\results\opto_sound_2025\context\noise\opto'];


[mouse_corr_stats, trial_indices_right_all, trial_indices_left_all] = neural_paircorr_acrossmice([1:24], params.info.mouse_date, params.info.serverid, dff_st,sound_sig_cells, [], 63:92,stim_trials_context, ctrl_trials_context,  corr_params,save_data_directory,pooled_cell_types);

%make cdf plot
p_values = cdf_corr_across_contexts(mouse_corr_stats, params.plot_info, save_data_directory,[-1:0.05:1]);

%%

snr_index_matrix = unpack_mouse_corr_stats(mouse_corr_stats, 'mean');
plot_info.y_lims = [-.9, .9];
% Set labels for plots.
plot_info.plot_labels = {'SNR','SNR'}; % Alternative could be {'Left Sounds','Right Sounds'}
plot_info.behavioral_contexts = {'Active','Passive'}; %decide which contexts to plot
overlap_labels = {'Active', 'Passive','Both'}; %{'Active', 'Passive','Both'}; % {'Active', 'Passive','Both'}; %{'Active', 'Passive','Spont','Both'}; %
params.plot_info = plot_info;
params.info.chosen_mice = [1:25];

%save directory
save_dir = ['V:\Connie\results\opto_sound_2025\context\snr\sounds'];% '/spont_sig'];% '/spont_sig']; %[info.savepath '/mod/' mod_params.mod_type '/spont_sig']; % Set directory to save figures.

snr_index_matrix = unpack_mouse_corr_stats(mouse_corr_stats, 'mean');
[snr_index_stats] = generate_mod_index_plots(snr_index_matrix, celltypes_ids, params, save_dir);