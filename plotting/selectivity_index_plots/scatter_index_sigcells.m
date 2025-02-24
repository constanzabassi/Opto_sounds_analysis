function modl_fit = scatter_index_sigcells(sig_mod_boot, all_celltypes, index, plot_info, save_path, string1, string2)
    % comes from scatter_passive_mod_opto_cells -
    % Plots scatter data for modulation indices of optogenetically modulated cells
    %
    % This function plots modulation index data for different cell types during passive
    % and optogenetic conditions. If `chosen_cells` is provided, only selected cells are used;
    % otherwise, all available cell data is plotted.
    %
    % Inputs:
    %   sig_mod_boot  - Cell array containing indices of selected cells for each mouse.
    %   all_celltypes - Cell array of structures with field names corresponding to different cell types.
    %   index         - Cell array containing modulation index data (two columns for x and y values).
    %   plot_info     - Structure containing color information for different cell types.
    %   save_path     - String specifying the directory for saving plots.
    %   string1       - Label for the x-axis.
    %   string2       - Label for the y-axis.
    %
    % Output:
    %   modl_fit      - Linear models fitted to the scatter plots for each cell type.
    
    % Define plot limits
    minmax = [-0.5, 0.5];

    % Create figure and set up the identity line
    figure(2324); clf;
    plot([minmax(1), minmax(2)], [minmax(1), minmax(2)], 'k--', 'Color', [0.5 0.5 0.5]);
    hold on;

    % Check if specific cells were chosen
    if ~isempty(sig_mod_boot)
        % Get field names for cell types (assumes consistent structure across mice)
        celltype_fields = fields(all_celltypes{1,1});

        % Loop over each cell type (e.g., pyramidal, SOM, PV)
        for cell_type = 1:3
            all_indices = []; % Collect indices for regression
            
            for dataset_index = 1:length(sig_mod_boot) % Iterate over mice
                % Get the indices of selected cells that match the current cell type
                selected_cells = sig_mod_boot{1, dataset_index}(ismember(sig_mod_boot{1, dataset_index}, all_celltypes{1, dataset_index}.(celltype_fields{cell_type})));

                % Scatter plot for the selected cells
                scatter(index{dataset_index,1}(selected_cells), index{dataset_index,2}(selected_cells), ...
                    'MarkerEdgeColor', plot_info.colors_celltypes(cell_type, :), ...
                    'MarkerEdgeAlpha', 1, 'LineWidth', 1.5);

                % Store indices for regression
                all_indices = [all_indices, [index{dataset_index,1}(selected_cells); index{dataset_index,2}(selected_cells)]];
            end

            % Fit a linear model and store the result
            modl_fit{cell_type} = fitlm(all_indices(1,:), all_indices(2,:));

            % Display R² value on the plot
            text(minmax(2), minmax(2) + 0.1 - 0.1 * cell_type, ...
                ['R² = ' num2str(modl_fit{cell_type}.Rsquared.Ordinary)], ...
                'Color', plot_info.colors_celltypes(cell_type, :), 'FontSize', 14);
        end
    else
        % If no specific cells were chosen, plot all available cell data
        celltype_fields = fields(all_celltypes{1,1});
        
        for cell_type = 1:3
            for dataset_index = 1:length(all_celltypes) % Iterate over mice
                scatter(index{dataset_index,1}(all_celltypes{1, dataset_index}.(celltype_fields{cell_type})), ...
                        index{dataset_index,2}(all_celltypes{1, dataset_index}.(celltype_fields{cell_type})), ...
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

    % Save figure if save path is provided
    if ~isempty(save_path)
        mkdir(save_path);
        cd(save_path);
        saveas(gcf, ['scatter_index_sigcells_' string1 '_' string2 '.png']);
        saveas(gcf, ['scatter_index_sigcells_' string1 '_' string2 '.svg']);
        exportgraphics(gcf,fullfile(save_dir, ['scatter_index_sigcells_' string1 '_' string2 '.pdf']), 'ContentType', 'vector');
    end
end
