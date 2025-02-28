function plot_side_preference(selectivity_results_all, params)
    pool_types = {'left', 'right'};
    
    figure('Position', [100 100 700 200]);
    
    for p_idx = 1:length(pool_types)
        pool = pool_types{p_idx};
        data = selectivity_results_all.both.(pool);
        
        % Create confusion matrix
        sides = {'left', 'right'};
        conf_mat = zeros(2);
        for i = 1:length(data.active_preferred)
            active_idx = strcmp(data.active_preferred{i}, sides);
            passive_idx = strcmp(data.passive_preferred{i}, sides);
            conf_mat(active_idx, passive_idx) = conf_mat(active_idx, passive_idx) + 1;
        end
        
        % Plot
        subplot(1, 2, p_idx);
        imagesc(conf_mat);
        colorbar;
        title([upper(pool) ' Selective']);
        xlabel('Passive Preferred');
        ylabel('Active Preferred');
        set(gca, 'XTick', 1:2, 'XTickLabel', sides);
        set(gca, 'YTick', 1:2, 'YTickLabel', sides);
        caxis([1,150])
        utils.set_current_fig;
        colorList= (colormaps.slanCM('plasma',100));
        colormap(colorList) % redblue
    end
end