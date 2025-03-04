function generate_neural_heatmaps(neural_structure, stim_trials_context, ctrl_trials_context, chosen_mice, params, type)
 %     % Generate heatmaps from neural_structure using z-scored data
%     %
%     % Inputs:
%     %   neural_structure - Structure with .stim and .ctrl fields [trials x neurons x frames]
%     %   trial_info      - Structure with trial conditions
%     %   current_dataset - Index of current dataset
%     %   context        - Context index (1=active, 2=passive)
%     %   params         - Plotting parameters
%  
    % Load trial information (adjust paths as needed) -virmen trial info left turns/sound condition/is stim
load('V:\Connie\results\opto_sound_2025\context\sound_info\active_all_trial_info_sounds.mat');
passive_all_trial_info_sounds = load('V:\Connie\results\opto_sound_2025\context\sound_info\passive_all_trial_info_sounds.mat').all_trial_info_sounds;

% Store data across contexts
    data_by_context = cell(1,2);  % To store {mean_stim/ctrl_left/right} for each context
    
    nContexts = size(stim_trials_context{1,1},2);
    for context = 1:nContexts
        % Initialize arrays for concatenating across datasets
        mean_stim_left = [];
        mean_stim_right = [];
        mean_ctrl_left = [];
        mean_ctrl_right = [];

        for current_dataset = chosen_mice
            % Get condition labels from trial info 
            if context == 1
                if strcmpi(type,'sound')
                    current_conditions = [all_trial_info_sounds(current_dataset).opto.condition];
                    current_conditions_ctrl = [all_trial_info_sounds(current_dataset).ctrl.condition,all_trial_info_sounds(current_dataset).sound_only.condition];
                else
                    current_conditions = [all_trial_info_sounds(current_dataset).opto.condition];
                    current_conditions_ctrl = [all_trial_info_sounds(current_dataset).ctrl.condition];
                end
            elseif context == 2
                if strcmpi(type,'sound')
                    current_conditions = [passive_all_trial_info_sounds(current_dataset).opto.condition];
                    current_conditions_ctrl = [passive_all_trial_info_sounds(current_dataset).ctrl.condition,passive_all_trial_info_sounds(current_dataset).sound_only.condition];
                else
                    current_conditions = [passive_all_trial_info_sounds(current_dataset).opto.condition];
                    current_conditions_ctrl = [passive_all_trial_info_sounds(current_dataset).ctrl.condition];
                end
            else
                current_conditions = [];
                current_conditions_ctrl = [];
            end
    
            % Get z-scored data and trial indices
            stim_data = neural_structure{1,current_dataset}.z_stim;
            ctrl_data = neural_structure{1,current_dataset}.z_ctrl;
            stim_trials = stim_trials_context{1, current_dataset}{1, context};
            ctrl_trials = ctrl_trials_context{1, current_dataset}{1, context};
    
            % Find trial types
            [~, ~, ~, ~, ~, ~, left_stim_all, left_ctrl_all, right_stim_all, right_ctrl_all] = ...
                find_sound_trials_single(stim_trials, ctrl_trials, current_conditions, current_conditions_ctrl);
    
            % Calculate and concatenate means
            mean_stim_left = [mean_stim_left; squeeze(mean(stim_data(left_stim_all,:,:), 1))];
            mean_stim_right = [mean_stim_right; squeeze(mean(stim_data(right_stim_all,:,:), 1))];
            mean_ctrl_left = [mean_ctrl_left; squeeze(mean(ctrl_data(left_ctrl_all,:,:), 1))];
            mean_ctrl_right = [mean_ctrl_right; squeeze(mean(ctrl_data(right_ctrl_all,:,:), 1))];
        end
        
        % Store data for this context
        data_by_context{context} = struct('stim_left', mean_stim_left, ...
                                        'stim_right', mean_stim_right, ...
                                        'ctrl_left', mean_ctrl_left, ...
                                        'ctrl_right', mean_ctrl_right);
    end
    
    params = get_heatmap_params();

    % Plot each context with specified sorting
    for context = 1:nContexts
        figure('Position', [100 100 800 500]);
        
        % Get data for current context
        curr_data = data_by_context{context};
        
        % Get sorting index if specified in params
        if isfield(params, 'sort_reference')
            ref = params.sort_reference;
            ref_context = ref.context;
            ref_data = data_by_context{ref_context};
            
            % Get reference data for sorting
            switch ref.condition
                case 'stim_left'
                    sort_data = ref_data.stim_left;
                case 'stim_right'
                    sort_data = ref_data.stim_right;
                case 'ctrl_left'
                    sort_data = ref_data.ctrl_left;
                case 'ctrl_right'
                    sort_data = ref_data.ctrl_right;
            end
            [sorted_idx, ~] = sort_neurons(sort_data, params);
        else
            % Default sorting using current context data
            [sorted_idx, ~] = sort_neurons(curr_data.stim_left, params);
        end
        
        % Plot using consistent sorting
        condition_fields = fields(curr_data);
        for current_condition = 1:length(condition_fields)
            subplot(2,2,current_condition)
            % Parse condition name into parts
            parts = strsplit(condition_fields{current_condition}, '_');
            if strcmpi(parts{1},'stim')
                title_str = sprintf('%s Sound with %s', upper(parts{2}), upper(parts{1}));
            else
                title_str = sprintf('%s Sound', (parts{2}));
            end

            plot_single_direction_heatmap(curr_data.(condition_fields{current_condition}), sorted_idx, title_str, params);
            utils.set_current_fig;
        end

        % Add common colorbar and save
        cb = colorbar;
        cb.Position = [0.92 0.1 0.02 0.8];
        cb.Label.String = 'Z-score';


    end

            %MAKE PLOTS SEPARATING STIM AND CONTROL
    % Create separate figures for stim and ctrl
    condition_types = {'stim', 'ctrl'};
    for cond = 1:length(condition_types)
        figure('Position', [100 100 600 500]);%figure('Position', [100 100 800 500]);
        current_cond = condition_types{cond};
        
        % Plot all contexts for left and right
        for direction = 1:2
            dir_name = {'Left', 'Right'};
            for context = 1:nContexts
                subplot(2, 2, (direction-1)*2 + context)
                field_name = sprintf('%s_%s', current_cond, lower(dir_name{direction}));

                [sorted_idx, ~] = sort_neurons(data_by_context{context}.(field_name), params)
                plot_single_direction_heatmap(data_by_context{context}.(field_name), ...
                    sorted_idx, sprintf('%s - %s Sound', ...
                    params.context_labels{context}, (dir_name{direction})), params);
                utils.set_current_fig;
                set(gca,'FontSize',12);

            end
        end
        
        % Add super title and colorbar
        sgtitle(sprintf('%s Responses', upper(current_cond)), 'FontWeight', 'normal');
        cb = colorbar;
        cb.Position = [0.9 0.3 0.02 0.3];%[0.92 0.1 0.02 0.8]; %left,bottom,idth,height
        cb.Label.String = 'Z-score';
        
        if isfield(params, 'savepath')
            saveas(gcf, fullfile(params.savepath, ...
                sprintf('%s_responses_heatmap.png', current_cond)));
            exportgraphics(gcf,fullfile(params.savepath, ...
                sprintf('%s_responses_heatmap.pdf', current_cond)), 'ContentType', 'vector');
        end
    end
        
        if isfield(params, 'savepath')
            saveas(gcf, fullfile(params.savepath, sprintf('context_%d_heatmap.png', context)));
        end
    end
% function generate_neural_heatmaps(neural_structure, stim_trials_context, ctrl_trials_context,chosen_mice, params, type)
%     % Generate heatmaps from neural_structure using z-scored data
%     %
%     % Inputs:
%     %   neural_structure - Structure with .stim and .ctrl fields [trials x neurons x frames]
%     %   trial_info      - Structure with trial conditions
%     %   current_dataset - Index of current dataset
%     %   context        - Context index (1=active, 2=passive)
%     %   params         - Plotting parameters
%  
%     % Load trial information (adjust paths as needed) -virmen trial info left turns/sound condition/is stim
% load('V:\Connie\results\opto_sound_2025\context\sound_info\active_all_trial_info_sounds.mat');
% passive_all_trial_info_sounds = load('V:\Connie\results\opto_sound_2025\context\sound_info\passive_all_trial_info_sounds.mat').all_trial_info_sounds;
% 
%     nContexts = size(stim_trials_context{1,1},2);
%     for context = 1:nContexts
%         % Initialize arrays for concatenating across datasets
%         mean_stim_left = [];
%         mean_stim_right = [];
%         mean_ctrl_left = [];
%         mean_ctrl_right = [];
% 
%         for current_dataset = chosen_mice
%             % Get condition labels from trial info 
%             if context == 1
%                 if strcmpi(type,'sound')%for sound alignment I included [control, sound_only trials] so I need to concatenate trial types here)
%                     current_conditions = [all_trial_info_sounds(current_dataset).opto.condition];
%                     current_conditions_ctrl = [all_trial_info_sounds(current_dataset).ctrl.condition,all_trial_info_sounds(current_dataset).sound_only.condition];
%                 else
%                     current_conditions = [all_trial_info_sounds(current_dataset).opto.condition];
%                     current_conditions_ctrl = [all_trial_info_sounds(current_dataset).ctrl.condition];
%                 end
%             elseif context == 2
%                 if strcmpi(type,'sound')%for sound alignment I included [control, sound_only trials] so I need to concatenate trial types here)
%                     current_conditions = [passive_all_trial_info_sounds(current_dataset).opto.condition];
%                     current_conditions_ctrl = [passive_all_trial_info_sounds(current_dataset).ctrl.condition,passive_all_trial_info_sounds(current_dataset).sound_only.condition];
%                 else
%                     current_conditions = [passive_all_trial_info_sounds(current_dataset).opto.condition];
%                     current_conditions_ctrl = [passive_all_trial_info_sounds(current_dataset).ctrl.condition];
%                 end
%             else %spont has no conditions
%                 current_conditions = [];
%                 current_conditions_ctrl = [];
%             end
%     
%             % Get z-scored stim and ctrl data
%             stim_data = neural_structure{1,current_dataset}.z_stim;
%             ctrl_data = neural_structure{1,current_dataset}.z_ctrl;
% 
%             % Get trial indices
%             stim_trials = stim_trials_context{1, current_dataset}{1, context};
%             ctrl_trials = ctrl_trials_context{1, current_dataset}{1, context};
%     
%             % Find left and right trials
%             [~, ~, ~, ~, ~, ~, left_stim_all, left_ctrl_all, right_stim_all, right_ctrl_all] = ...
%                 find_sound_trials_single(stim_trials, ctrl_trials, current_conditions, current_conditions_ctrl);
%     
%             % Calculate mean responses for this dataset
%             curr_mean_stim_left = squeeze(mean(stim_data(left_stim_all,:,:), 1));
%             curr_mean_stim_right = squeeze(mean(stim_data(right_stim_all,:,:), 1));
%             curr_mean_ctrl_left = squeeze(mean(ctrl_data(left_ctrl_all,:,:), 1));
%             curr_mean_ctrl_right = squeeze(mean(ctrl_data(right_ctrl_all,:,:), 1));
%             
%             % Concatenate across datasets
%             mean_stim_left = [mean_stim_left; curr_mean_stim_left];
%             mean_stim_right = [mean_stim_right; curr_mean_stim_right];
%             mean_ctrl_left = [mean_ctrl_left; curr_mean_ctrl_left];
%             mean_ctrl_right = [mean_ctrl_right; curr_mean_ctrl_right];
%         end
% 
%         % Create figure for this context
%         figure('Position', [100 100 1200 800]);
%         
%         % Get sorting index based on combined responses
%         if ~strcmpi(type, 'sound')
%             combined_responses = [mean_stim_left, mean_stim_right];
%         else
%             combined_responses = [mean_ctrl_left, mean_ctrl_right];
%         end
%         params = get_heatmap_params();
%         [sorted_idx, ~] = sort_neurons(combined_responses, params);
%         
%         % Plot stimulation trials
%         subplot(2,2,1)
%         [sorted_idx, ~] = sort_neurons(mean_stim_left, params);
%         plot_single_direction_heatmap(mean_stim_left, sorted_idx, 'Left Sound - Stim', params);
%         utils.set_current_fig;
%         subplot(2,2,2)
%         [sorted_idx, ~] = sort_neurons(mean_stim_right, params);
%         plot_single_direction_heatmap(mean_stim_right, sorted_idx, 'Right Sound - Stim', params);
%         utils.set_current_fig;
%         
%         % Plot control trials
%         subplot(2,2,3)
%         [sorted_idx, ~] = sort_neurons(mean_ctrl_left, params);
%         plot_single_direction_heatmap(mean_ctrl_left, sorted_idx, 'Left Sound - Ctrl', params);
%         utils.set_current_fig;
%         subplot(2,2,4)
%         [sorted_idx, ~] = sort_neurons(mean_ctrl_right, params);
%         plot_single_direction_heatmap(mean_ctrl_right, sorted_idx, 'Right Sound - Ctrl', params);
%         utils.set_current_fig;
% %         sgtitle(sprintf('%s Context - All Datasets (n=%d neurons)', ...
% %             params.context_labels{context}, size(mean_stim_left,1)), ...
% %             'FontWeight', 'normal');
%         
%         % Add common colorbar
%         cb = colorbar;
%         cb.Position = [0.92 0.1 0.02 0.8];
%         cb.Label.String = 'Z-score';
%         
%         if isfield(params, 'savepath')
%             saveas(gcf, fullfile(params.savepath, ...
%                 sprintf('context_%d_heatmap.png', context)));
%         end
%     end
% end
