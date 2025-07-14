%% load context data

load('V:\Connie\results\opto_sound_2025\context\data_info\all_celltypes.mat');
load('V:\Connie\results\opto_sound_2025\context\data_info\context_data.mat');
plot_info = plotting_config(); %plotting params

keep context_data all_celltypes
%% define axis
[proj,proj_ctrl,proj_norm,proj_norm_ctrl, weights,trial_corr_context,percent_correct,act,act_norm_ctrl,act_norm,percent_correct_concat,proj_concat,proj_concat_norm] = find_axis(context_data.dff, [1:24], all_celltypes,[],[],[]); %,{50:59,63:73}

save_dir = 'V:\Connie\results\opto_sound_2025\context\axis_lme_plots\dff';

%% plot mean projection traces across datasets
celltype = 4; %4 = all

plot_proj_mean_traces([1:24],proj_ctrl, 'sound',celltype, [61:62],[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);
plot_proj_mean_traces([1:24],proj_ctrl, 'context',celltype, [61:62],[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);
plot_proj_mean_traces([1:24],proj, 'stim',celltype, [61:62],[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);

frames_to_avg = 50:59;
bin_edges = [-2:0.4:2];%
hist_stats =  histogram_axis_across_contexts([1:24],proj_ctrl, 'context',celltype, bin_edges,frames_to_avg,[0,0,0;.5,.5,.5],{'Active','Passive'},save_dir);

%% model
celltype = 4;
frame_range_pre= 50:59;
frame_range_post = 63:93;

%sounds(predicted) vs context
[lme_sound,tbl_sound,proj_all_sound,engagement_proj_all_sound,context_all_sound] = mixed_linear_model(proj_ctrl, 'Sound',celltype,frame_range_pre,frame_range_post);
[lme_sound_pass,tbl_sound_pass,~,~,~] = mixed_linear_model(proj_ctrl, 'Sound',celltype,frame_range_pre,frame_range_post,'Context',1);

%stim(predicted) vs context
[lme_stim,tbl_stim,proj_all_stim,engagement_proj_all_stim,context_all_stim] = mixed_linear_model(proj, 'Stim',celltype,frame_range_pre,frame_range_post);
[lme_stim_pass,tbl_stim_pass,~,~,~] = mixed_linear_model(proj, 'Stim',celltype,frame_range_pre,frame_range_post,'Context',1);

% %stim(predicted) vs sound
% [lme_stim_sound,tbl_stim_sound,proj_all_stim_sound,engagement_proj_all_stim_sound,context_all_stim_sound] = mixed_linear_model(proj, 'Stim',celltype,frame_range_post,frame_range_post,'Sound');
% [lme_stim_sound_pass,tbl_stim_sound_pass,~,~,~] = mixed_linear_model(proj, 'Stim',celltype,frame_range_post,frame_range_post,'Sound',1);

%sound(predicted) vs stim
[lme_sound_stim,tbl_sound_stim,proj_all_sound_stim,engagement_proj_all_sound_stim,context_all_sound_stim] = mixed_linear_model(proj,'Sound' ,celltype,frame_range_post,frame_range_post,'Stim');
[lme_sound_stim_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Sound' ,celltype,frame_range_post,frame_range_post,'Stim',1);


%% make plots
coeffs_stim = extract_and_rename_coefficients(lme_stim, lme_stim_pass, 'passive', ...
    {'Engagement Slope (A)', 'Engagement Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});

coeffs_sound = extract_and_rename_coefficients(lme_sound, lme_sound_pass, 'passive', ...
    {'Engagement Slope (A)', 'Engagement Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});

% coeffs_stim_sound = extract_and_rename_coefficients(lme_stim_sound, lme_stim_sound_pass, 'passive', ...
%     {'Sound Slope (A)', 'Sound Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});

coeffs_sound_stim = extract_and_rename_coefficients(lme_sound_stim, lme_sound_stim_pass, 'passive', ...
    {'Stim Slope (A)', 'Stim Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});


plot_fixed_effects(coeffs_sound,coeffs_stim,  save_dir, [0.3,0.2,0.6 ; 1,0.7,0],[]); %“Active Engagement Effect”“Passive Engagement Effect”“Context Offset (Passive - Active)”“Interaction (Slope Difference)”%{'Intercept','Engagement effect(A)','P vs A offset', 'Context:Engagement (P vs A)'}; %{"Effect of engagement (active)", "Passive vs. active shift", "Change in engagement effect (passive vs. active)"}
% plot_fixed_effects(coeffs_stim_sound, coeffs_stim_sound, save_dir, [0,0,0],[]);
plot_fixed_effects(coeffs_sound_stim, coeffs_sound_stim, save_dir, [0,0,0],[]);

% plot_fixed_effects(lme_stim, lme_sound, save_dir, [0.3,0.2,0.6 ; 1,0.7,0],[]); %“Active Engagement Effect”“Passive Engagement Effect”“Context Offset (Passive - Active)”“Interaction (Slope Difference)”%{'Intercept','Engagement effect(A)','P vs A offset', 'Context:Engagement (P vs A)'}; %{"Effect of engagement (active)", "Passive vs. active shift", "Change in engagement effect (passive vs. active)"}
% plot_fixed_effects(lme_stim_sound, lme_stim_sound, save_dir, [0,0,0],[]);


p_val_sound = plot_me_residuals(tbl_sound,'Sound',save_dir);
p_val_stim = plot_me_residuals(tbl_stim,'Stim',save_dir);

plot_me_regression_lines(lme_sound,engagement_proj_all_sound,proj_all_sound,context_all_sound,'Sound Projection',save_dir);
plot_me_regression_lines(lme_stim,engagement_proj_all_stim,proj_all_stim,context_all_stim,'Stim Projection',save_dir);

% plot_me_regression_lines(lme_stim_sound,engagement_proj_all_stim_sound,proj_all_stim_sound,context_all_stim_sound,'Stim Projection',save_dir,'Sound');
plot_me_regression_lines(lme_sound_stim,engagement_proj_all_sound_stim,proj_all_sound_stim,context_all_sound_stim,'Sound Projection',save_dir,'Stim');
%% performance
lme_perf = plot_errorbar_performance_lme(percent_correct_concat,engagement_proj_all_sound,engagement_proj_all_stim,context_all_sound,context_all_stim,save_dir);

%% weights
colors_medium = [0.37 0.75 0.49 %green
                0.17 0.35 0.8  %blue
                0.82 0.04 0.04];
edges_values_weights = [-.1,.1];
num_bins_weights = 20;
[weight_all_celltype,weight_ct_stats] = histogram_weights_celltypes_vs_axis([1:24],weights, 'Context' ,all_celltypes, edges_values_weights,num_bins_weights,colors_medium,save_dir);

%% extra plots
plot_me_violin_random_intercept(proj_all_sound,context_all_sound,lme_sound, 'Sound');
plot_me_violin_random_intercept(proj_all_stim,context_all_stim,lme_stim, 'Stim');


plot_error_bars_response_vs_axis_allbinned(proj,proj_all_stim,engagement_proj_all_stim,context_all_stim,'Stim',[0,0,0;0.5,0.5,0.5],[]);
plot_error_bars_response_vs_axis_allbinned(proj_ctrl,proj_all_sound,engagement_proj_all_sound,context_all_sound,'Sound',[0,0,0;0.5,0.5,0.5],[]);
%% cell type effects (what is being predicted) vs same projections as before (using all cell types)
frame_range_pre= 50:59;
frame_range_post = 63:93;
for celltype = 1:3
    %update directory
    save_dir_celltype = strcat(save_dir,upper(plot_info.celltype_names{celltype}));

    %sounds(predicted) vs context
    [lme_sound,tbl_sound,proj_all_sound,engagement_proj_all_sound,context_all_sound] = mixed_linear_model(proj_ctrl, 'Sound',[celltype,4],frame_range_pre,frame_range_post);
    [lme_sound_pass,tbl_sound_pass,~,~,~] = mixed_linear_model(proj_ctrl, 'Sound',[celltype,4],frame_range_pre,frame_range_post,'Context',1);
    
    %stim(predicted) vs context
    [lme_stim,tbl_stim,proj_all_stim,engagement_proj_all_stim,context_all_stim] = mixed_linear_model(proj, 'Stim',[celltype,4],frame_range_pre,frame_range_post);
    [lme_stim_pass,tbl_stim_pass,~,~,~] = mixed_linear_model(proj, 'Stim',[celltype,4],frame_range_pre,frame_range_post,'Context',1);
        
    %sound(predicted) vs stim
    [lme_sound_stim,tbl_sound_stim,proj_all_sound_stim,engagement_proj_all_sound_stim,context_all_sound_stim] = mixed_linear_model(proj,'Sound' ,[celltype,4],frame_range_post,frame_range_post,'Stim');
    [lme_sound_stim_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Sound' ,[celltype,4],frame_range_post,frame_range_post,'Stim',1);

    coeffs_stim = extract_and_rename_coefficients(lme_stim, lme_stim_pass, 'passive', ...
    {'Engagement Slope (A)', 'Engagement Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});

    coeffs_sound = extract_and_rename_coefficients(lme_sound, lme_sound_pass, 'passive', ...
        {'Engagement Slope (A)', 'Engagement Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});
    
    coeffs_sound_stim = extract_and_rename_coefficients(lme_sound_stim, lme_sound_stim_pass, 'passive', ...
        {'Stim Slope (A)', 'Stim Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});

    %main plots
    plot_fixed_effects(coeffs_sound,coeffs_stim,  save_dir_celltype, [0.3,0.2,0.6 ; 1,0.7,0],[]); %“Active Engagement Effect”“Passive Engagement Effect”“Context Offset (Passive - Active)”“Interaction (Slope Difference)”%{'Intercept','Engagement effect(A)','P vs A offset', 'Context:Engagement (P vs A)'}; %{"Effect of engagement (active)", "Passive vs. active shift", "Change in engagement effect (passive vs. active)"}
    plot_fixed_effects(coeffs_sound_stim, coeffs_sound_stim, save_dir_celltype, [0,0,0],[]);
    p_val_sound = plot_me_residuals(tbl_sound,'Sound',save_dir_celltype);
    p_val_stim = plot_me_residuals(tbl_stim,'Stim',save_dir_celltype);
    
    plot_me_regression_lines(lme_sound,engagement_proj_all_sound,proj_all_sound,context_all_sound,'Sound Projection',save_dir_celltype);
    plot_me_regression_lines(lme_stim,engagement_proj_all_stim,proj_all_stim,context_all_stim,'Stim Projection',save_dir_celltype);
    
    plot_me_regression_lines(lme_sound_stim,engagement_proj_all_sound_stim,proj_all_sound_stim,context_all_sound_stim,'Sound Projection',save_dir_celltype,'Stim');
end
%% look at inhibitory neurons (decide which one to use)
%sound(predicted-PYR) vs inhibitory neurons (using same axis just pre vs post)
[lme_sound_pv,tbl_sound_stim,proj_all_sound_pv,engagement_proj_all_sound_pv,context_all_sound_pv] = mixed_linear_model(proj,'Sound' ,[1,3],frame_range_pre,frame_range_post,'Sound');
[lme_sound_pv_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Sound' ,[1,3],frame_range_pre,frame_range_post,'Sound',1);
[lme_sound_som,tbl_sound_stim,proj_all_sound_som,engagement_proj_all_sound_som,context_all_sound_som] = mixed_linear_model(proj,'Sound' ,[1,2],frame_range_pre,frame_range_post,'Sound');
[lme_sound_som_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Sound' ,[1,2],frame_range_pre,frame_range_post,'Sound',1);
%stim(predicted-PYR) vs inhibitory neurons (using same axis just pre vs post)
[lme_sound_pv,tbl_sound_stim,proj_all_sound_pv,engagement_proj_all_sound_pv,context_all_sound_pv] = mixed_linear_model(proj,'Stim' ,[1,3],frame_range_pre,frame_range_post,'Stim');
[lme_sound_pv_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Stim' ,[1,3],frame_range_pre,frame_range_post,'Stim',1);
[lme_sound_som,tbl_sound_stim,proj_all_sound_som,engagement_proj_all_sound_som,context_all_sound_som] = mixed_linear_model(proj,'Stim' ,[1,2],frame_range_pre,frame_range_post,'Stim');
[lme_sound_som_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Stim' ,[1,2],frame_range_pre,frame_range_post,'Stim',1);

%orthogonalized IH projections 
[lme_sound_pv,tbl_sound_stim,proj_all_sound_pv,engagement_proj_all_sound_pv,context_all_sound_pv] = mixed_linear_model(proj,'Stim' ,[1,3],frame_range_post,frame_range_post,'Stim_post');
[lme_sound_pv_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Stim' ,[1,3],frame_range_post,frame_range_post,'Stim_post',1);
[lme_sound_som,tbl_sound_stim,proj_all_sound_som,engagement_proj_all_sound_som,context_all_sound_som] = mixed_linear_model(proj,'Stim' ,[1,2],frame_range_post,frame_range_post,'Stim_post');
[lme_sound_som_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Stim' ,[1,2],frame_range_post,frame_range_post,'Stim_post',1);

[lme_sound_pv,tbl_sound_stim,proj_all_sound_pv,engagement_proj_all_sound_pv,context_all_sound_pv] = mixed_linear_model(proj,'Stim' ,[1,3],frame_range_pre,frame_range_post,'Stim_pre');
[lme_sound_pv_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Stim' ,[1,3],frame_range_pre,frame_range_post,'Stim_pre',1);
[lme_sound_som,tbl_sound_stim,proj_all_sound_som,engagement_proj_all_sound_som,context_all_sound_som] = mixed_linear_model(proj,'Stim' ,[1,2],frame_range_pre,frame_range_post,'Stim_pre');
[lme_sound_som_pass,tbl_sound_stim_pass,~,~,~] = mixed_linear_model(proj,'Stim' ,[1,2],frame_range_pre,frame_range_post,'Stim_pre',1);


% organize outputs
coeffs_sound_ih = extract_and_rename_coefficients(lme_sound_pv, lme_sound_pv_pass, 'passive', ...
    {'IH Slope (A)', 'IH Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});
coeffs_sound_som = extract_and_rename_coefficients(lme_sound_som, lme_sound_som_pass, 'passive', ...
    {'IH Slope (A)', 'IH Slope (P)', 'Context Diff. (A - P)', 'Slope Diff. (A - P)'});


%inhibitory effects
plot_fixed_effects(coeffs_sound_som, coeffs_sound_ih, save_dir, plot_info.colors_celltypes(2:3,:),[],'sound_vs_preIH');
plot_me_regression_lines(lme_sound_pv,engagement_proj_all_sound_pv,proj_all_sound_pv,context_all_sound_pv,'Sound Projection',save_dir,'PV');
plot_me_regression_lines(lme_sound_som,engagement_proj_all_sound_som,proj_all_sound_som,context_all_sound_som,'Sound Projection',save_dir,'SOM');

plot_fixed_effects(coeffs_sound_som, coeffs_sound_ih, save_dir, plot_info.colors_celltypes(2:3,:),[],'stim_vs_preIH');
plot_me_regression_lines(lme_sound_pv,engagement_proj_all_sound_pv,proj_all_sound_pv,context_all_sound_pv,'Stim Projection',save_dir,'PV');
plot_me_regression_lines(lme_sound_som,engagement_proj_all_sound_som,proj_all_sound_som,context_all_sound_som,'Stim Projection',save_dir,'SOM');

plot_fixed_effects(coeffs_sound_som, coeffs_sound_ih, save_dir, plot_info.colors_celltypes(2:3,:),[],'stim_vs_postIH');
plot_me_regression_lines(lme_sound_pv,engagement_proj_all_sound_pv,proj_all_sound_pv,context_all_sound_pv,'Stim Projection',save_dir,'PV (post)');