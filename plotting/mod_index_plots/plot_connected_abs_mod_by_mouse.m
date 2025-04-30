function mod_stats = plot_connected_abs_mod_by_mouse(save_dir, mod_index_by_dataset, mouseID, plot_info, varargin)
    % Plot modulation index with connected lines for each mouse
    % Inputs:
    %   save_dir - directory to save results
    %   mod_index_by_dataset - cell array {dataset, context, celltype}
    %   mouseID - array indicating mouse ID for each dataset
    %   plot_info - structure with plotting parameters
    
    figure(700);clf
    
    num_contexts = 2;
    unique_mice = unique(mouseID);
    n_mice = length(unique_mice);
    n_celltypes =  size( mod_index_by_dataset,3);;
    count = 0;
    x_lines = 0:num_contexts*n_celltypes+1;
    
    if nargin > 4
        ylims = varargin{1,1};
    end
    
    for celltype = 1:n_celltypes
        mean_cell_all = [];
        
        % Get x positions for this celltype
        x_pos = x_lines((celltype-1)*num_contexts+2:celltype*num_contexts+1);
        
        % Calculate and plot individual mouse lines
        for m = 1:n_mice
            curr_mouse = unique_mice(m);
            mouse_datasets = find(mouseID == curr_mouse);
            
            % Get mean for each context
            mouse_means = zeros(1, num_contexts);
            for context = 1:num_contexts
                % Average across datasets for this mouse
                mouse_data = [];
                for d = mouse_datasets
                    curr_data = abs(mod_index_by_dataset{d,context,celltype});
                    mouse_data = [mouse_data; curr_data(:)];
                end
                mouse_means(context) = nanmean(mouse_data);
            end
            
            % Plot connected line for this mouse
            plot([x_pos(1)+.2,x_pos(2)-.2], mouse_means, '-', 'Color', ...
                [plot_info.colors_celltypes(celltype,:), 0.3], 'LineWidth', 1)
            hold on
            
            % Store mouse means for statistics
            for context = 1:num_contexts
                mod_stats.stats(celltype,context).mouse_means(m) = mouse_means(context);
            end
        end
        
        % Calculate and plot mean ± SEM across mice for each context
        for context = 1:num_contexts
            count = count + 1;
            curr_means = mod_stats.stats(celltype,context).mouse_means;
%             mean_cel = nanmean(curr_means);
%             err = std(curr_means, 'omitnan') / sqrt(sum(~isnan(curr_means)));

            valid_data = curr_means(~isnan(curr_means));
            if ~isempty(valid_data)
                mean_cel = mean(valid_data);
                err = std(valid_data) / sqrt(length(valid_data));
                
                % Store stats with number of valid datasets
                mod_stats.stats(celltype,context).valid_means = valid_data;
                mod_stats.stats(celltype,context).mean = mean_cel;
                mod_stats.stats(celltype,context).sem = err;
                mod_stats.stats(celltype,context).n_valid_datasets = length(valid_data);
                 % Calculate 95% CI using bootstrapping
                mod_stats.stats(celltype,context).basic_stats =  get_basic_stats(valid_data)
            end

%             % Plot error bar (CI)
%             errorbar(x_pos(context), mean_cel, err_low, err_up, 'o', ...
%                 'Color', plot_info.colors_celltypes(celltype,:), ...
%                 'LineWidth', 1.3, 'MarkerSize', 3,'MarkerFaceColor', plot_info.colors_celltypes(celltype,:))
            
            % Plot error bar (SEM)
            errorbar(x_pos(context), mean_cel, err, 'o', ...
                'Color', plot_info.colors_celltypes(celltype,:), ...
                'LineWidth', 1, 'MarkerSize', 2,'MarkerFaceColor', plot_info.colors_celltypes(celltype,:))
            
            mean_cell_all = [mean_cell_all, mean_cel];
        end
        
        % Statistical testing for this cell type
        ct = 0;
        possible_tests = nchoosek(1:num_contexts,2);
        if max(mean_cell_all) > 0.1
            y_val = max(mean_cell_all) + 0.03;
        else
            y_val = max(mean_cell_all);
        end
        
        for t = 1:size(possible_tests,1)
            data1 = mod_stats.stats(celltype,possible_tests(t,1)).valid_means;
            data2 = mod_stats.stats(celltype,possible_tests(t,2)).valid_means;
            
            [p_val_mod(t,celltype), ~, effectsize(t,celltype)] = permutationTest_updatedcb(...
                data1, data2, 10000, 'paired', 1);
            
            if p_val_mod(t,celltype) < 0.05/n_celltypes
                xline_vars = possible_tests(t,:) + ((celltype-1)*num_contexts);
                ct = ct + 0.03;
                plot_pval_star(0, y_val+ct, p_val_mod(t,celltype), xline_vars, ...
                    0.01, plot_info.colors_celltypes(celltype,:))
            end
        end
    end
    
    % Format plot
    xlim([x_lines(1) x_lines(end)])
    xticks(x_lines(2:end-1))
    xticklabels(repmat(plot_info.behavioral_contexts, 1, n_celltypes))
    ylabel({'Absolute Modulation';'Index'})
    
    % Set axis limits
    if nargin > 4
        ylim(ylims);
    else
        yli = ylim;
        ylim([0,yli(2)]);
    end
    
    % Set figure properties
    set(gca,'FontSize',12);
    set(gcf,'Color','w')
    set(gca,'FontName','Arial')
%     set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
    box off
    xtickangle(45)
    set(gcf,'units','points','position',[10,100,(400/n_celltypes*num_contexts),170]);
    utils.set_current_fig;
    
    % Store statistics
    mod_stats.tests = possible_tests;
    mod_stats.p_test = 'paired permutation across mice';
    mod_stats.p_val_mod = p_val_mod;
    mod_stats.effectsize = effectsize;
    mod_stats.n_mice = n_mice;
    
    % Save results
    if ~isempty(save_dir)
        mkdir(save_dir)
        cd(save_dir)
        saveas(700,strcat('abs_mod_index_connected_lines_n',num2str(n_mice),'.svg'));
        saveas(700,strcat('abs_mod_index_connected_lines_n',num2str(n_mice),'.fig'));
        exportgraphics(figure(700),strcat('abs_mod_index_connected_lines_n',num2str(n_mice),'_datasets.pdf'), 'ContentType', 'vector');
        save(strcat('abs_mod_index_stats_connected_lines_n',num2str(n_mice)),'mod_stats');
    end
end