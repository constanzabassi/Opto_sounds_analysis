load('V:\Connie\results\opto_sound_2025\context\data_info\all_celltypes.mat');
load('V:\Connie\results\opto_sound_2025\context\data_info\context_data.mat');
plot_info = plotting_config(); %plotting params

keep context_data all_celltypes plot_info
%% define axis
% [proj,proj_ctrl,proj_norm,proj_norm_ctrl, weights,trial_corr_context,percent_correct,act,act_norm_ctrl,act_norm,percent_correct_concat,proj_concat,proj_concat_norm] = find_axis_updated(context_data.dff, [1:24], all_celltypes,[]); %,{50:59,63:73}

split_params.divisions = 3; split_params.random_or_not = 0; split_params.splits = 3;
[proj,proj_ctrl,proj_norm,proj_norm_ctrl, weights,trial_corr_context,percent_correct,act,act_norm_ctrl,act_norm,percent_correct_concat,proj_concat,proj_concat_norm] = find_axis_updated_specify_splits(context_data.dff, [1:24], all_celltypes,[],split_params); %,{50:59,63:73}

save_dir = 'W:\Connie\results\Bassi2025\fig4\updated_3cv\';%'V:\Connie\results\opto_sound_2025\context\axis_lme_plots_updated\dff';
save_dir2 = 'W:\Connie\results\Bassi2025\fig5\updated_3cv\';%'V:\Connie\results\opto_sound_2025\context\axis_lme_plots_updated\dff';

%% plot mean projection traces across datasets (finds means across splits first
celltype = 4; %4 = all
plot_proj_meansplits_traces([1:24],proj_ctrl, 'sound',celltype, [61:62],[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);
plot_proj_meansplits_traces([1:24],proj_ctrl, 'context',celltype, [61:62],[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);
plot_proj_meansplits_traces([1:24],proj, 'stim',celltype, [61:62],[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);

% plot_proj_mean_traces([1:24],squeeze(proj_ctrl(1,:,:,:)), 'sound',celltype, [61:62],[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);

frames_to_avg = 50:59;
bin_edges = [-2:0.4:2];%
hist_stats =  histogram_axis_across_contexts_splits([1:24],proj_ctrl, 'context',celltype, bin_edges,frames_to_avg,[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);

frames_to_avg = 63:93;
bin_edges = [-2:0.4:2];%
hist_stats =  histogram_axis_across_contexts_splits([1:24],proj_ctrl, 'Sound',celltype, bin_edges,frames_to_avg,[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);

frames_to_avg = 63:93;
bin_edges = [-1.5:0.4:2.5];%
hist_stats =  histogram_axis_across_contexts_splits([1:24],proj, 'Stim',celltype, bin_edges,frames_to_avg,[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);
%% model
celltype = 4;
frame_range_pre= 50:59;
frame_range_post = 63:93;
%sound (predicted) vs engagement axis
[lme_sound,tbl_sound,proj_all_sound,engagement_proj_all_sound,context_all_sound,corr_mean, corr_all] = ...
    linear_regression_corr_model(proj_norm_ctrl, 'Sound',celltype,frame_range_pre,frame_range_post,[1:2]);
%stim(predicted) vs engagement axis
[lme_stim,tbl_stim,proj_all_stim,engagement_proj_all_stim,context_all_stim,corr_mean_stim, corr_all_stim] = ...
    linear_regression_corr_model(proj_norm, 'Stim',celltype,frame_range_pre,frame_range_post,[1:2]);

% % %sound(predicted) vs stim
[lme_sound_stim,tbl_sound_stim,proj_all_sound_stim,engagement_proj_all_sound_stim,context_all_sound_stim,corr_mean_sound_stim, corr_all_sound_stim] = ...
    linear_regression_corr_model(proj_norm,'Sound' ,celltype,frame_range_post,frame_range_post,[1],'Stim');

[lme_sound_stim_pass,tbl_sound_stim_pass,~,~,context_all_sound_stim_pass,corr_mean_sound_stim_pass, corr_all_sound_stim_pass] = linear_regression_corr_model(proj_norm,'Sound' ,celltype,frame_range_post,frame_range_post,[2],'Stim');

%% scatter plots of trials and linear regression lines
plot_linear_regression_lines(lme_sound,tbl_sound,context_all_sound,'Sound Projection',save_dir);
plot_linear_regression_lines(lme_stim,tbl_stim,context_all_stim,'Stim Projection',save_dir);
plot_linear_regression_lines(lme_sound_stim,tbl_sound_stim,context_all_sound_stim,'Sound Projection',save_dir2,'Stim');
plot_linear_regression_lines(lme_sound_stim_pass,tbl_sound_stim_pass,context_all_sound_stim_pass,'Sound Projection',save_dir2,'Stim');

%% performance vs engagment!
nDatasets = 24;
for d = 1:nDatasets
    percent_correct_all{d} = horzcat(percent_correct_concat{:,d});   % concatenate across splits
    %find mean across neurons
    concat_activity{d} = vertcat(act{:,d,4}); %trials x neurons
    mean_act = mean(concat_activity{d}(:,50:59),2);
    activity_all{d} = mean_act';
end
plot_performance_vs_engagement_axis(percent_correct_all,engagement_proj_all_sound,engagement_proj_all_stim,context_all_sound,context_all_stim,[10,5],save_dir,[-2,2]);
plot_performance_vs_engagement_axis(activity_all,engagement_proj_all_sound,engagement_proj_all_stim,context_all_sound,context_all_stim,[10,5],save_dir,[-2,2]);
%% plot weights across cell types
colors_medium = [0.37 0.75 0.49 %green
                0.17 0.35 0.8  %blue
                0.82 0.04 0.04];
edges_values_weights = [-.1,.1];
num_bins_weights = 20;
[weight_all_celltype,weight_ct_stats] = histogram_weights_celltypes_vs_axis_splits([1:24],weights, 'Context' ,all_celltypes, edges_values_weights,num_bins_weights,colors_medium,save_dir);

%% using random splits (choose first one to show)
% model correlations across splits
celltype = 4;
frame_range_pre= 50:59;
frame_range_post = 63:93;
%sound (predicted) vs engagement axis
[~,~,~,~,~,corr_mean, corr_all, corr_stats] = ...
    linear_regression_corr_model(proj_norm_ctrl, 'Sound',celltype,frame_range_pre,frame_range_post,[1:2]);
%stim(predicted) vs engagement axis
[~,~,~,~,~,corr_mean_stim, corr_all_stim, corr_stats_stim] = ...
    linear_regression_corr_model(proj_norm, 'Stim',celltype,frame_range_pre,frame_range_post,[1:2]);

% % %sound(predicted) vs stim
[~,~,~,~,~,corr_mean_sound_stim, corr_all_sound_stim, corr_stats_sound_stim] = ...
    linear_regression_corr_model(proj_norm,'Sound' ,celltype,frame_range_post,frame_range_post,[1],'Stim');

[~,~,~,~,~,corr_mean_sound_stim_pass, corr_all_sound_stim_pass, corr_stats_sound_stim_pass] = linear_regression_corr_model(proj_norm,'Sound' ,celltype,frame_range_post,frame_range_post,[2],'Stim');

%% from random have to choose example split
celltype = 4;
frame_range_pre= 50:59;
frame_range_post = 63:93;
example_split = 1;

%sound (predicted) vs engagement axis
[lme_sound,tbl_sound,proj_all_sound,engagement_proj_all_sound,context_all_sound] = ...
    linear_regression_corr_model(proj_norm_ctrl(example_split,:,:,:), 'Sound',celltype,frame_range_pre,frame_range_post,[1:2]);
%stim(predicted) vs engagement axis
[lme_stim,tbl_stim,proj_all_stim,engagement_proj_all_stim,context_all_stim,corr_mean_stim] = ...
    linear_regression_corr_model(proj_norm(example_split,:,:,:), 'Stim',celltype,frame_range_pre,frame_range_post,[1:2]);

% % %sound(predicted) vs stim
[lme_sound_stim,tbl_sound_stim,proj_all_sound_stim,engagement_proj_all_sound_stim,context_all_sound_stim] = ...
    linear_regression_corr_model(proj_norm(example_split,:,:,:),'Sound' ,celltype,frame_range_post,frame_range_post,[1],'Stim');

[lme_sound_stim_pass,tbl_sound_stim_pass,~,~,context_all_sound_stim_pass,corr_mean_sound_stim_pass] = linear_regression_corr_model(proj_norm(example_split,:,:,:),'Sound' ,celltype,frame_range_post,frame_range_post,[2],'Stim');

% scatter plots of trials and linear regression lines for example session
plot_linear_regression_lines(lme_sound,tbl_sound,context_all_sound,'Sound Projection',save_dir);
plot_linear_regression_lines(lme_stim,tbl_stim,context_all_stim,'Stim Projection',save_dir);
plot_linear_regression_lines(lme_sound_stim,tbl_sound_stim,context_all_sound_stim,'Sound Projection',save_dir2,'Stim');
plot_linear_regression_lines(lme_sound_stim_pass,tbl_sound_stim_pass,context_all_sound_stim_pass,'Sound Projection',save_dir2,'Stim');

colors_medium = [0.37 0.75 0.49 %green
                0.17 0.35 0.8  %blue
                0.82 0.04 0.04];
edges_values_weights = [-.1,.1];
num_bins_weights = 20;
[weight_all_celltype,weight_ct_stats] = histogram_weights_celltypes_vs_axis_splits([1:24],weights, 'Context' ,all_celltypes, edges_values_weights,num_bins_weights,colors_medium,save_dir);
