function plot_selectivity_comparison(selectivity_results_all, savepath)
    % Plot active vs passive modulation for each selectivity pool
    
    % Setup
    pool_types = {'left', 'right', 'nonsel'};
    directions = {'left_mod', 'right_mod'};
    plot_titles = {'Left Sound Modulation', 'Right Sound Modulation'};
    contexts = {'active', 'passive'};
    
    % Create figure for each pool type
    for p_idx = 1:length(pool_types)
        pool = pool_types{p_idx};
        
        % Create figure
        fig = figure('Position', [100 100 900 400]);
        
        % Create subplots for left and right modulation
        for d_idx = 1:length(directions)
            direction = directions{d_idx};
            subplot(1, 2, d_idx);
            
            % Get data for active and passive contexts from 'both' selective cells
            active_data = selectivity_results_all.both.(pool).(['active_' direction]);
            passive_data = selectivity_results_all.both.(pool).(['passive_' direction]);
            
            % Create matrix for heatmap [neurons x contexts]
            heatmap_data = [active_data, passive_data];
            
            % Sort based on average response
            if p_idx == 1
                sorting_ids = sort_avg_context_index(heatmap_data);
            elseif p_idx == 2
                active_data = selectivity_results_all.both.right.(['active_right_mod']);
                passive_data = selectivity_results_all.both.right.(['passive_right_mod']);
                sorting_ids = sort_avg_context_index([active_data, passive_data]);
            else
                sorting_ids = sort_avg_context_index(heatmap_data);
            end
            
            % Create heatmap
            make_heatmap_sorted(heatmap_data, ...
                contexts, ...                    % x-label (contexts)
                'Neurons', ...                   % y-label
                [-0.4,0.4], ...     % color limits
                sorting_ids);                    % sorting indices
            
            % Format plot
            title(plot_titles{d_idx});
            colorbar;
            
            % Add statistics if needed
            [~, p] = ttest2(active_data, passive_data);
            if p < 0.05
                text(0.5, -0.1, sprintf('p = %.3f', p), ...
                    'Units', 'normalized', 'HorizontalAlignment', 'center');
            end

        end
        
        % Add overall title
        sgtitle(sprintf('%s Selective Cells (n=%d)', ...
            upper(pool), length(selectivity_results_all.both.(pool).cell_indices)), ...
            'FontSize', 14);
        
        % Save figure if path provided
        if ~isempty(savepath)
            mkdir(savepath)
            saveas(fig, fullfile(savepath, ...
                sprintf('modulation_comparison_%s_selective.png', pool)));
            saveas(fig, fullfile(savepath, ...
                sprintf('modulation_comparison_%s_selective.fig', pool)));
        end
    end
end
