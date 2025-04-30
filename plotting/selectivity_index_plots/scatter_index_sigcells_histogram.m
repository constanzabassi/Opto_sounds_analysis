function modl_fit = scatter_index_sigcells_histogram(sig_mod_boot, all_celltypes, index, plot_info, save_path, string1, string2, varargin)
    if nargin > 7
        minmax = varargin{1,1};
    else
        minmax = [-1, 1];
    end

    % Create figure
    figure(2324); clf;
    
    % Define main axes
    mainAx = axes('Position', [0.15 0.15 0.65 0.65]);
    hold(mainAx, 'on');
    plot(mainAx, [minmax(1), minmax(2)], [minmax(1), minmax(2)], '--', 'Color', [0.5 0.5 0.5]);

    % Axes for histograms
    topAx = axes('Position', [0.15 0.81 0.65 0.15]); % X-axis histogram
    rightAx = axes('Position', [0.81 0.15 0.15 0.65]); % Y-axis histogram

    hold(topAx, 'on'); hold(rightAx, 'on');
    axis(topAx, 'tight'); axis(rightAx, 'tight');
    topAx.XTick = []; rightAx.YTick = [];
    set(topAx, 'XLim', minmax, 'YTick', []);
    set(rightAx, 'YLim', minmax, 'XTick', []);

    modl_fit = cell(1, 3);
    celltype_fields = fields(all_celltypes{1});
    n_celltypes = length(celltype_fields);

    for cell_type = 1:n_celltypes
        all_x = []; all_y = [];
        for dataset_index = 1:length(all_celltypes)
            if ~isempty(sig_mod_boot)
                selected_cells = sig_mod_boot{dataset_index}(ismember(sig_mod_boot{dataset_index}, ...
                    all_celltypes{dataset_index}.(celltype_fields{cell_type})));
            else
                selected_cells = all_celltypes{dataset_index}.(celltype_fields{cell_type});
            end
            x = index{dataset_index,1}(selected_cells);
            y = index{dataset_index,2}(selected_cells);

            % Scatter plot
            scatter(mainAx, x, y, 'MarkerEdgeColor', plot_info.colors_celltypes(cell_type, :), ...
                'MarkerEdgeAlpha', .7, 'LineWidth', 1.5);

            % Collect data for fitting and histograms
            all_x = [all_x; x(:)];
            all_y = [all_y; y(:)];
        end

        if ~isempty(all_x)
            modl_fit{cell_type} = fitlm(all_x, all_y);
            text(mainAx, minmax(1)+0.05, minmax(2) +0.05 - 0.1 * cell_type, ...
                sprintf('R² = %.3f', modl_fit{cell_type}.Rsquared.Ordinary), ...
                'Color', plot_info.colors_celltypes(cell_type, :), 'FontSize', 14);

            % Plot histograms
            histogram(topAx, all_x, 'BinLimits', minmax, 'FaceColor', plot_info.colors_celltypes(cell_type, :), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Orientation', 'vertical');

            histogram(rightAx, all_y, 'BinLimits', minmax, 'FaceColor', plot_info.colors_celltypes(cell_type, :), ...
                'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Orientation', 'horizontal');
        end
    end


    % Finalize plot
    xlabel(mainAx, string1);
    ylabel(mainAx, string2);
    xlim(mainAx, minmax);
    ylim(mainAx, minmax);
    box(mainAx, 'off');
    utils.set_current_fig;


    % Save figure
    sig_cel_string = ~isempty(sig_mod_boot);
    if ~isempty(save_path)
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        saveas(gcf, fullfile(save_path, ['scatter_index_sigcells' num2str(sig_cel_string) '_' string1 '_' string2 '.png']));
        saveas(gcf, fullfile(save_path, ['scatter_index_sigcells' num2str(sig_cel_string) '_' string1 '_' string2 '.svg']));
        exportgraphics(gcf, fullfile(save_path, ['scatter_index_sigcells' num2str(sig_cel_string) '_' string1 '_' string2 '.pdf']), 'ContentType', 'vector');
    end
end
