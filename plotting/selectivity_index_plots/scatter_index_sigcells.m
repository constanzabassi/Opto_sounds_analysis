function modl_fit = scatter_index_sigcells(sig_mod_boot, all_celltypes, index, plot_info, save_path, string1, string2)
    % Plots scatter data for modulation indices of optogenetically modulated cells
    % 
    % This function visualizes modulation index data for different neuron types under passive 
    % and optogenetic conditions. If `sig_mod_boot` is provided, only selected cells are plotted.
    %
    % Inputs:
    %   sig_mod_boot  - Cell array of indices for selected neurons (each entry corresponds to a dataset).
    %   all_celltypes - Cell array of structures containing different neuron types.
    %   index         - Cell array with modulation index data (two columns: x-values and y-values).
    %   plot_info     - Struct containing colors for different cell types.
    %   save_path     - Directory path for saving figures.
    %   string1       - X-axis label.
    %   string2       - Y-axis label.
    %
    % Output:
    %   modl_fit      - Cell array containing fitted linear models for each cell type.

    % Define plot limits
    minmax = [-0.5, 0.5];

    % Create figure and plot identity line
    figure(2324); clf;
    plot([minmax(1), minmax(2)], [minmax(1), minmax(2)], '--', 'Color', [0.5 0.5 0.5]);
    hold on;

    % Initialize model fit results
    modl_fit = cell(1, 3);  

    % Determine if we're plotting selected or all cells
    if ~isempty(sig_mod_boot)
        celltype_fields = fields(all_celltypes{1,1}); % Get neuron type labels

        for cell_type = 1:3
            all_indices = []; % Store all data points for regression
            
            for dataset_index = 1:length(sig_mod_boot) % Iterate over datasets
                % Identify selected cells matching the current neuron type
                selected_cells = sig_mod_boot{dataset_index}(ismember(sig_mod_boot{dataset_index}, all_celltypes{dataset_index}.(celltype_fields{cell_type})));

                % Scatter plot of selected cells
                scatter(index{dataset_index,1}(selected_cells), index{dataset_index,2}(selected_cells), ...
                    'MarkerEdgeColor', plot_info.colors_celltypes(cell_type, :), ...
                    'MarkerEdgeAlpha', 1, 'LineWidth', 1.5);

                % Store data points for regression
                all_indices = [all_indices; index{dataset_index,1}(selected_cells), index{dataset_index,2}(selected_cells)];
            end

            % Fit a linear model and store results
            if ~isempty(all_indices)
                modl_fit{cell_type} = fitlm(all_indices(:,1), all_indices(:,2));

                % Display R² value
                text(minmax(2), minmax(2) + 0.1 - 0.1 * cell_type, ...
                    sprintf('R² = %.3f', modl_fit{cell_type}.Rsquared.Ordinary), ...
                    'Color', plot_info.colors_celltypes(cell_type, :), 'FontSize', 14);
            end
        end
    else
        % Plot all available neuron data if no specific cells were chosen
        celltype_fields = fields(all_celltypes{1,1});
        
        for cell_type = 1:3
            for dataset_index = 1:length(all_celltypes)
                scatter(index{dataset_index,1}(all_celltypes{dataset_index}.(celltype_fields{cell_type})), ...
                        index{dataset_index,2}(all_celltypes{dataset_index}.(celltype_fields{cell_type})), ...
                        'MarkerEdgeColor', plot_info.colors_celltypes(cell_type, :), ...
                        'MarkerEdgeAlpha', 1, 'LineWidth', 1.5);
            end
        end
    end

    % Finalize plot settings
    hold off;
    box off;
    xlim(minmax);
    ylim(minmax);
    xlabel(string1);
    ylabel(string2);

    % Format figure
    utils.set_current_fig;

    % Save figure if a save path is provided
    if ~isempty(save_path)
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        
        saveas(gcf, fullfile(save_path, ['scatter_index_sigcells_' string1 '_' string2 '.png']));
        saveas(gcf, fullfile(save_path, ['scatter_index_sigcells_' string1 '_' string2 '.svg']));
        exportgraphics(gcf, fullfile(save_path, ['scatter_index_sigcells_' string1 '_' string2 '.pdf']), 'ContentType', 'vector');
    end
end
