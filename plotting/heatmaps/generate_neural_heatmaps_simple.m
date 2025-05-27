function generate_neural_heatmaps_simple(neural_structure, stim_trials_context, ctrl_trials_context,sig_mod_boot, chosen_mice, params,type, context_to_plot,varargin)
 %     % Generate heatmaps from neural_structure using z-scored data
%     %
%     % Inputs:
%     %   neural_structure - Structure with .stim and .ctrl fields [trials x neurons x frames]
%     %   trial_info      - Structure with trial conditions
%     %   current_dataset - Index of current dataset
%     %   context        - Context index (1=active, 2=passive)
%     %   params         - Plotting parameters
%  
nContexts = length(context_to_plot); %size(stim_trials_context{1,1},2);

% Store data across contexts
data_by_context = cell(1,nContexts);  % To store for each context
data_by_context_ctrl = cell(1,nContexts);  % To store 

params = get_heatmap_params();
positions = utils.calculateFigurePositions(1, 6, .5, []);

if nargin > 8
    minmax = varargin{1,1};
else
    minmax = [-.5,1];
end

% Plot each context with specified sorting
context_count = 0;
for context = context_to_plot
    context_count = context_count+1;
    
    temp = [];temp2 = [];
    for dataset_index = chosen_mice
        cellCount = size(neural_structure{1,dataset_index}.z_stim,2);
        if ~isempty(sig_mod_boot)
            mod_cells = sig_mod_boot{1,dataset_index}; %
            save_name = 'sig_cells';
        else
            mod_cells = 1:cellCount;
            save_name = 'all_cells';
        end
        temp = [temp;squeeze(mean(neural_structure{1,dataset_index}.z_stim(stim_trials_context{1,dataset_index}{1,context},mod_cells,:)))];
        temp2 = [temp2;squeeze(mean(neural_structure{1,dataset_index}.z_ctrl(ctrl_trials_context{1,dataset_index}{1,context},mod_cells,:)))];

    end
    data_by_context{context_count} = temp;
    data_by_context_ctrl{context_count} = temp2;


    %sort by first context
    if context_count == 1
        % Get sorting index if specified in params
        sort_data = temp;
        sort_data2 = temp2;
        
        % Get reference data for sorting
        [sorted_idx, ~] = sort_neurons(sort_data, params);
        [sorted_idx2, ~] = sort_neurons(sort_data2, params);
    end
    
end

figure(1);clf;
%PLOT BASIC HEATMAPS
for context = 1:nContexts
    % Plot using consistent sorting
    subplot(1,nContexts,context)

    plot_single_direction_heatmap(data_by_context{context}, sorted_idx, [params.context_labels{context_to_plot(context)}], params);
    utils.set_current_fig;

    % Add common colorbar and save
    if context == nContexts
    cb = colorbar;
%     cb.Position = [0.92 0.1 0.02 0.8];
    cb.Label.String = 'Z-scored ΔF/F';
    caxis(minmax);
    end
    set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(context, :));
end
if isfield(params, 'savepath')
    saveas(gcf, fullfile(params.savepath, ...
        sprintf('heatmap_stim_%s_%s_%d.png',type,save_name,context_to_plot)));
    exportgraphics(gcf,fullfile(params.savepath, ...
        sprintf('heatmap_stim_%s_%s_%d.pdf', type,save_name,context_to_plot)), 'ContentType', 'vector');
end

figure(2);clf;
%PLOT BASIC HEATMAPS
for context = 1:nContexts
    % Plot using consistent sorting
    subplot(1,nContexts,context)

    plot_single_direction_heatmap(data_by_context_ctrl{context}, sorted_idx2, [params.context_labels{context_to_plot(context)}], params);
    utils.set_current_fig;

    % Add common colorbar and save
    if context == nContexts
    cb = colorbar;
%     cb.Position = [0.92 0.1 0.02 0.8];
    cb.Label.String = 'Z-scored ΔF/F';
    caxis(minmax);
    end
    set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(context, :));
end
if isfield(params, 'savepath')
    saveas(gcf, fullfile(params.savepath, ...
        sprintf('heatmap_ctrl_%s_%s_%d.png',type,save_name,context_to_plot)));
    exportgraphics(gcf,fullfile(params.savepath, ...
        sprintf('heatmap_ctrl_%s_%s_%d.pdf', type,save_name,context_to_plot)), 'ContentType', 'vector');
end

end
