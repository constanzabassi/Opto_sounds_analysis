function [proj,proj_ctrl,proj_norm,proj_ctrl_norm, weights,trial_corr_context,percent_correct,real_activity_all,real_activity_all_ctrl,percent_correct_concat,proj_concat] = find_axis(dff_st, chosen_mice, all_celltypes,sig_mod_boot,sig_mod_boot2,sig_mod_boot3,varargin)
        total_trials = {};
        possible_celltypes = fieldnames(all_celltypes{1,1});

        %LOAD VIRMEN TRIAL INFO
        all_trial_info = load('V:\Connie\results\opto_sound_2025\context\sound_info\active_all_trial_info_sounds.mat').all_trial_info_sounds; %all_trial_info

        rng(2025);
    % Loop through the selected mice datasets
    for current_dataset = chosen_mice
        current_dataset
        % Example: Randomly split trials into train (80%) and test (20%)
        total_trials_stim = [size(dff_st{1, current_dataset}.stim,1),size(dff_st{2, current_dataset}.stim,1)]; %get total # trials across contexts
        total_trials_ctrl = [size(dff_st{1, current_dataset}.ctrl,1),size(dff_st{2, current_dataset}.ctrl,1)]; %get total # trials across contexts

        %separate trials threeway
        ctx_count_stim = 0;
        ctx_count_ctrl = 0;
        for ctx = 1:2
            % Stim trials
            [trainA_stim{ctx}, trainB_stim{ctx}, test_stim{ctx}] = ...
                split_trials_threeway(total_trials_stim(ctx));
            
            % Add context offset
            trainA_stim{ctx} = trainA_stim{ctx} + ctx_count_stim;
            trainB_stim{ctx} = trainB_stim{ctx} + ctx_count_stim;
            test_stim{ctx}   = test_stim{ctx}   + ctx_count_stim;
            
            % Ctrl trials
            [trainA_ctrl{ctx}, trainB_ctrl{ctx}, test_ctrl{ctx}] = ...
                split_trials_threeway(total_trials_ctrl(ctx));
            
            % Add context offset
            trainA_ctrl{ctx} = trainA_ctrl{ctx} + ctx_count_ctrl;
            trainB_ctrl{ctx} = trainB_ctrl{ctx} + ctx_count_ctrl;
            test_ctrl{ctx}   = test_ctrl{ctx}   + ctx_count_ctrl;
            
            ctx_count_stim = ctx_count_stim + total_trials_stim(ctx);
            ctx_count_ctrl = ctx_count_ctrl + total_trials_ctrl(ctx);
        end
        
        % Combine across contexts if needed
        trainA_stim_all = [trainA_stim{1}, trainA_stim{2}];
        trainB_stim_all = [trainB_stim{1}, trainB_stim{2}];
        test_stim_all   = [test_stim{1}, test_stim{2}];

        trainA_ctrl_all = [trainA_ctrl{1}, trainA_ctrl{2}];
        trainB_ctrl_all = [trainB_ctrl{1}, trainB_ctrl{2}];
        test_ctrl_all   = [test_ctrl{1}, test_ctrl{2}];


        for celltype = 1:4

            if celltype == 4
                mod_cells = 1:size(dff_st{1, current_dataset}.stim,2);
                if ~isempty(sig_mod_boot)
                    mod_cells = sig_mod_boot{current_dataset};
                end
            else
                mod_cells = all_celltypes{1,current_dataset}.(possible_celltypes{celltype});
                if ~isempty(sig_mod_boot)
                    mod_cells = mod_cells(find(ismember(mod_cells,sig_mod_boot{ current_dataset})));
                end
            end
            % Get the data for the stim and control contexts
            stim_matrix = [dff_st{1, current_dataset}.stim(:,mod_cells,:); dff_st{2, current_dataset}.stim(:,mod_cells,:)]; % stim + sound (active and passive)
            ctrl_matrix = [dff_st{1, current_dataset}.ctrl(:,mod_cells,:); dff_st{2, current_dataset}.ctrl(:,mod_cells,:)]; % sound only (active and passive)
            [~, n_neurons, frames] = size(stim_matrix);
%             [trials_ctrl] = size(ctrl_matrix, 1);

            % Define the frames for after and before stimulus
            aframes = 63:93;  % After stimulus
            bframes = 50:59;   % Before stimulus
    
            % --- 1) Calculate Sound Axis (pre-post sound) ---
            data_matrix = ctrl_matrix;  % Use the control matrix for the sound axis calculation
            prepost_trial = [];
            tr_ct = 0;
            for tr = trainA_ctrl_all;%train_trials_ctr %1:trials_ctrl
                tr_ct = tr_ct +1;
                prepost_trial(tr,:,:) = [nanmean(data_matrix(tr, :, aframes), [3]);
                                       nanmean(data_matrix(tr, :, bframes), [3])]; % [post, pre]
            end
            sound_diff = mean((squeeze(prepost_trial(:, 1,:))) -squeeze(prepost_trial(:, 2,:)),[1]);%mean((squeeze(prepost_trial(:, 1,:))));% -squeeze(prepost_trial(:, 2,:)),[1]); %difference size is neurons
%             sound_diff = mean((squeeze(prepost_trial(:, 1,:))));% just taking post response
            norm_sound_diff = sound_diff ./ sqrt(sum(sound_diff.^2)); % Normalize sound axis
    

            if celltype == 4
                mod_cells2 = 1:size(dff_st{1, current_dataset}.stim,2);
                if ~isempty(sig_mod_boot2)
                    mod_cells2 = sig_mod_boot2{current_dataset};
                end
            else
                mod_cells2 = all_celltypes{1,current_dataset}.(possible_celltypes{celltype});
                if ~isempty(sig_mod_boot2)
                    mod_cells2 = mod_cells2(find(ismember(mod_cells2,sig_mod_boot2{current_dataset})));
                    if isempty(mod_cells2)
                        continue
                    end
                end
            end

            stim_matrix2 = [dff_st{1, current_dataset}.stim(:,mod_cells2,:); dff_st{2, current_dataset}.stim(:,mod_cells2,:)]; % stim + sound (active and passive)
            ctrl_matrix2 = [dff_st{1, current_dataset}.ctrl(:,mod_cells2,:); dff_st{2, current_dataset}.ctrl(:,mod_cells2,:)]; % sound only (active and passive)
            
            % --- 2) Calculate Stimulus + Sound Axis (stim + sound - sound) ---
            stim_sound_mean = mean(stim_matrix2(trainA_stim_all,:, aframes), [1, 3]);
            sound_mean = mean(ctrl_matrix2(trainA_ctrl_all, :, aframes), [1, 3]);
            stim_axis = stim_sound_mean - sound_mean; % Difference between stim + sound and sound
            norm_stim_diff = stim_axis ./ sqrt(sum(stim_axis.^2)); % Normalize stimulus axis
  
    
            % --- 3) Calculate Active-Passive Axis (active - passive) ---

            %use different trials to find this axis (just in case)
            if celltype == 4
                mod_cells3 = 1:size(dff_st{1, current_dataset}.stim,2);
                if ~isempty(sig_mod_boot3)
                    mod_cells3 = sig_mod_boot3{ current_dataset};
                end
            else
                mod_cells3 = all_celltypes{1,current_dataset}.(possible_celltypes{celltype});
                if ~isempty(sig_mod_boot3)
                    mod_cells3 = mod_cells3(find(ismember(mod_cells3,sig_mod_boot3{current_dataset})));
                end
            end
            stim_matrix3 = [dff_st{1, current_dataset}.stim(:,mod_cells3,:); dff_st{2, current_dataset}.stim(:,mod_cells3,:)]; % stim + sound (active and passive)
            ctrl_matrix3 = [dff_st{1, current_dataset}.ctrl(:,mod_cells3,:); dff_st{2, current_dataset}.ctrl(:,mod_cells3,:)]; % sound only (active and passive)

            active_passive_sound_mean = [
            nanmean(ctrl_matrix3(trainB_ctrl_all(find(ismember( trainB_ctrl_all,trainB_ctrl{1}))), :, bframes), [1, 3]); % active sound total_trials{current_dataset, 1, 2}
            nanmean(ctrl_matrix3(trainB_ctrl_all(find(ismember( trainB_ctrl_all,trainB_ctrl{2}))), :, bframes), [1, 3])]; % passive sound total_trials{current_dataset, 2, 2}
            context_axis = active_passive_sound_mean(1,:) - active_passive_sound_mean(2,:); % active - passive
            norm_context_diff = context_axis ./ sqrt(sum(context_axis.^2)); % Normalize context axis
            % do pre of stim trials also (to increase number of trials for
            % behavioral performance plots)
            active_passive_sound_mean_stim = [
            nanmean(stim_matrix3(trainB_stim_all(find(ismember( trainB_stim_all,trainB_stim{1}))), :, bframes), [1, 3]); % active sound total_trials{current_dataset, 1, 2}
            nanmean(stim_matrix3(trainB_stim_all(find(ismember( trainB_stim_all,trainB_stim{2}))), :, bframes), [1, 3])]; % passive sound total_trials{current_dataset, 2, 2}
            context_axis_stim = active_passive_sound_mean_stim(1,:) - active_passive_sound_mean_stim(2,:); % active - passive
            norm_context_diff_stim = context_axis_stim ./ sqrt(sum(context_axis_stim.^2)); % Normalize context axis

    
            %5) PROJECT THE DATA! onto test trials!!
            %also save real activity in case we want to look at it
            
            tr = 0;
            for trial = test_ctrl_all 
                tr = tr+1;
                sound_proj_ctrl(tr,:) = squeeze(ctrl_matrix(trial,:,:))'*norm_sound_diff';
                stim_proj_ctrl(tr,:) = squeeze(ctrl_matrix2(trial,:,:))'*norm_stim_diff';
%                 context_proj_ctrl(tr,:) = squeeze(ctrl_matrix(trial,:,:))'*norm_context_diff';
                real_activity_ctrl(tr,:) = squeeze(mean(ctrl_matrix3(trial,:,:),2))'; %find mean across cells
            end

            tr = 0;
            for trial = test_stim_all%1: trials_stim
                tr = tr+1;
                sound_proj_stim(tr,:) = squeeze(stim_matrix(trial,:,:))'*norm_sound_diff';
                stim_proj_stim(tr,:) = squeeze(stim_matrix2(trial,:,:))'*norm_stim_diff';
%                 context_proj_stim(tr,:) = squeeze(stim_matrix(trial,:,:))'*norm_context_diff';
                real_activity_stim(tr,:) = squeeze(mean(stim_matrix3(trial,:,:),2))'; %find mean across cells
            end

            tr = 0;
            for trial = test_ctrl_all 
                tr = tr+1;
                context_proj_ctrl(tr,:) = squeeze(ctrl_matrix3(trial,:,:))'*norm_context_diff';
            end

            tr = 0;
            for trial = test_stim_all
                tr = tr+1;
                context_proj_stim(tr,:) = squeeze(stim_matrix3(trial,:,:))'*norm_context_diff_stim';
            end

            %axis save into a structure?
            weights{current_dataset,celltype}.sound = norm_sound_diff;
            weights{current_dataset,celltype}.stim = norm_stim_diff;
            weights{current_dataset,celltype}.context = norm_context_diff;

            %projections
            %stim
            proj_concat{current_dataset,celltype,1}.sound = zscore(sound_proj_ctrl);
            proj_concat{current_dataset,celltype,1}.stim = zscore(stim_proj_stim);
            proj_concat{current_dataset,celltype,1}.context_stim = zscore(context_proj_stim);
            proj_concat{current_dataset,celltype,1}.context_sound = zscore(context_proj_ctrl);
            for context = 1:2
                proj{current_dataset,celltype,context}.sound = sound_proj_stim(find(ismember( test_stim_all,test_stim{context})),:); %total_trials{current_dataset, context, 1}
                proj{current_dataset,celltype,context}.stim = stim_proj_stim(find(ismember( test_stim_all,test_stim{context})),:);
%                 proj{current_dataset,celltype,context}.context = context_proj_stim(find(ismember( test_trials,test_stim_trials_idx{context})),:);
                proj{current_dataset,celltype,context}.context = context_proj_stim(find(ismember( test_stim_all,test_stim{context})),:);
               
                %ctrl
                proj_ctrl{current_dataset,celltype,context}.sound = sound_proj_ctrl(find(ismember( test_ctrl_all ,test_ctrl{context})),:);
                proj_ctrl{current_dataset,celltype,context}.stim = stim_proj_ctrl(find(ismember( test_ctrl_all ,test_ctrl{context})),:);
%                 proj_ctrl{current_dataset,celltype,context}.context = context_proj_ctrl(find(ismember( test_ctrl_all ,test_ctrl_trials_idx{context})),:);
                proj_ctrl{current_dataset,celltype,context}.context = context_proj_ctrl(find(ismember( test_ctrl_all ,test_ctrl{context})),:);

                %concat stim and ctrl for more trials
                proj_concat{current_dataset,celltype,context}.context = zscore([proj{current_dataset,celltype,context}.context;proj_ctrl{current_dataset,celltype,context}.context]);

                % find correlations before and after
                trial_corr_context{current_dataset,celltype,context}.sound =  corr(nanmean(sound_proj_ctrl(find(ismember( test_ctrl_all ,test_ctrl{context})),aframes),2),nanmean(context_proj_ctrl(find(ismember( test_ctrl_all ,test_ctrl{context})),bframes),2),'Type','Pearson');% corr(nanmean(sound_proj_ctrl(total_trials{current_dataset, context, 2},:),2),nanmean(context_proj_ctrl(total_trials{current_dataset, context, 2},:),2),'Type','Pearson');
                trial_corr_context{current_dataset,celltype,context}.stim = corr(nanmean(stim_proj_stim(find(ismember( test_stim_all,test_stim{context})),aframes),2),nanmean(context_proj_stim(find(ismember( test_stim_all,test_stim{context})),bframes),2),'Type','Pearson');

                 % Normalize (z-score) the projections per dataset/context
                for key = ["sound", "stim", "context"]
                    this_proj = proj_ctrl{current_dataset,celltype,context}.(key);
                    proj_ctrl_norm{current_dataset,celltype,context}.(key) = (this_proj - mean(this_proj,2)) ./ std(this_proj,[],2);
                
                    this_proj = proj{current_dataset,celltype,context}.(key);
                    proj_norm{current_dataset,celltype,context}.(key) = (this_proj - mean(this_proj,2)) ./ std(this_proj,[],2);
                end

                %also save real actity
                real_activity_all_ctrl{current_dataset,celltype,context}.sound = zscore(real_activity_ctrl(find(ismember( test_ctrl_all ,test_ctrl{context})),:),[],2);
                real_activity_all_ctrl{current_dataset,celltype,context}.context = zscore(real_activity_ctrl(find(ismember( test_ctrl_all ,test_ctrl{context})),:),[],2);
                real_activity_all{current_dataset,celltype,context}.stim = zscore(real_activity_stim(find(ismember( test_stim_all,test_stim{context})),:),[],2);
                real_activity_all{current_dataset,celltype,context}.context = zscore(real_activity_stim(find(ismember( test_stim_all,test_stim{context})),:),[],2);


            end

            %save correct trials
            all_correct_trials = [all_trial_info(current_dataset).ctrl(:).correct];
            percent_correct{current_dataset} = all_correct_trials(test_ctrl_all (find(ismember( test_ctrl_all ,test_ctrl{1}))));

            %concatenate ctrl and stim to have more trials (performance
            %only)
            all_correct_trials_stim = [all_trial_info(current_dataset).opto(:).correct];
            percent_correct_concat{current_dataset} = [all_correct_trials_stim(test_stim_all (find(ismember( test_stim_all ,test_stim{1})))),percent_correct{current_dataset}];
%     

        end
    end
end