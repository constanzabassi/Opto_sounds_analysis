function [mouse_corr_stats, trial_indices_right_all, trial_indices_left_all] = neural_paircorr_acrossmice(chosen_mice, mouse_date, server, dff_st,combined_thres, frames_before, frames_after, stim_trials_context, ctrl_trials_context, function_params, save_data_directory)
    % Initialize output variables
    mouse_corr_stats = struct(); %cell(length(chosen_mice), 4);
    trial_indices_right_all = cell(length(chosen_mice), 1);
    trial_indices_left_all = cell(length(chosen_mice), 1);
    
    % Load cell type IDs if needed
    load('V:\Connie\results\opto_2024\context\data_info\all_celltypes.mat');

    %load trial info
    load('V:\Connie\results\opto_sound_2025\context\sound_info\active_all_trial_info_sounds.mat');
    passive_all_trial_info_sounds = load('V:\Connie\results\opto_sound_2025\context\sound_info\passive_all_trial_info_sounds.mat').all_trial_info_sounds;


        %Process across contexts
    for context_idx = 1:2
        context = context_idx ;%function_params.context(context_idx);
    
        % Process each mouse
        for current_dataset_index = 1:length(chosen_mice)
            current_dataset_index

            % Get mouse info
            current_mouse = chosen_mice(current_dataset_index);
            current_mouse_date = mouse_date{current_dataset_index};
            current_server = server{current_dataset_index};
                            
            % Get trial indices based on context
            % Get condition labels from trial info 
            if context == 1
                if strcmpi(function_params.type,'sound')
                    current_conditions = [all_trial_info_sounds(current_dataset_index).opto.condition];
                    current_conditions_ctrl = [all_trial_info_sounds(current_dataset_index).ctrl.condition,all_trial_info_sounds(current_dataset_index).sound_only.condition];
                else
                    current_conditions = [all_trial_info_sounds(current_dataset_index).opto.condition];
                    current_conditions_ctrl = [all_trial_info_sounds(current_dataset_index).ctrl.condition];
                end
            elseif context == 2
                if strcmpi(function_params.type,'sound')
                    current_conditions = [passive_all_trial_info_sounds(current_dataset_index).opto.condition];
                    current_conditions_ctrl = [passive_all_trial_info_sounds(current_dataset_index).ctrl.condition,passive_all_trial_info_sounds(current_dataset_index).sound_only.condition];
                else
                    current_conditions = [passive_all_trial_info_sounds(current_dataset_index).opto.condition];
                    current_conditions_ctrl = [passive_all_trial_info_sounds(current_dataset_index).ctrl.condition];
                end
            else
                current_conditions = [];
                current_conditions_ctrl = [];
            end
    
            % Get z-scored data and trial indices
%             stim_data = dff_st{1,current_dataset_index}.z_stim;
%             ctrl_data = dff_st{1,current_dataset_index}.z_ctrl;
            stim_trials = stim_trials_context{1, current_dataset_index}{1, context};
            ctrl_trials = ctrl_trials_context{1, current_dataset_index}{1, context};
    
            % Find trial types
            [~, ~, ~, ~, ~, ~, left_stim_all, left_ctrl_all, right_stim_all, right_ctrl_all] = ...
                find_sound_trials_single(stim_trials, ctrl_trials, current_conditions, current_conditions_ctrl);
            
            % Store trial indices
            if strcmpi(function_params.type,'sound')
                trial_indices_left_all{current_dataset_index} = [left_ctrl_all];
                trials_left = left_ctrl_all;
                trial_indices_right_all{current_dataset_index} = [right_ctrl_all];
                trials_right = right_ctrl_all;
            else
                trial_indices_left_all{current_dataset_index} = [left_stim_all];
                trials_left = [left_stim_all];
                trial_indices_right_all{current_dataset_index} = [right_stim_all];
                trials_right = [right_stim_all];
            end

            % Get current dataset's cell types
            function_params.current_celltypes = all_celltypes{1, current_dataset_index};
            celtype_fields  = fields(all_celltypes{1, 1});

            if strcmp(function_params.calc_type,'corr')
                for celltype = 1:3 
                    celltype_field = celtype_fields{celltype};
                    % Process correlations
                    [corr_left, corr_right] = process_trials_neural_corr(...
                        chosen_mice, current_dataset_index, mouse_date, dff_st, ...
                        frames_before, frames_after, trials_left, trials_right, combined_thres, celltype_field , function_params, save_data_directory);
                    
    
                    mouse_corr_stats.(celltype_field){current_dataset_index, context_idx, 1} = corr_left;
                    mouse_corr_stats.(celltype_field){current_dataset_index, context_idx, 2} = corr_right;
                end
            else
                  [corr_left, corr_right] = process_trials_neural_snr(...
                        chosen_mice, current_dataset_index, mouse_date, dff_st, ...
                        frames_before, frames_after, trials_left, trials_right, combined_thres, function_params, save_data_directory);

                    mouse_corr_stats.snr{current_dataset_index, context_idx, 1} = corr_left;
                    mouse_corr_stats.snr{current_dataset_index, context_idx, 2} = corr_right;

            end

        end
    end
end
