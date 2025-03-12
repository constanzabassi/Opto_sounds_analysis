corr_params.corr_bins = 3; %how many frames to bin
corr_params.nonoverlap_bins = 1;
corr_params.type = 'stim';
corr_params.field_name_neural = 'stim'; %or sound
save_data_directory =['V:\Connie\results\opto_sound_2025\context\corr\opto'];


[mouse_corr_stats, trial_indices_right_all, trial_indices_left_all] = neural_paircorr_acrossmice([1:24], params.info.mouse_date, params.info.serverid, dff_st,combined_sig_cells, 50:59, 63:92,stim_trials_context, ctrl_trials_context,  corr_params,save_data_directory);

%make cdf plot
p_values = cdf_corr_across_contexts(mouse_corr_stats, params.plot_info, save_data_directory,[-.2:0.05:.2]);