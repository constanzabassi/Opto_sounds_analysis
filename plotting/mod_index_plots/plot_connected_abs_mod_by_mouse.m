function mod_stats = plot_connected_abs_mod_by_mouse(save_dir, mod_index_by_dataset, mouseID, plot_info, varargin)
    % Plot modulation index with connected lines for each mouse
    % Inputs:
    %   save_dir - directory to save results
    %   mod_index_by_dataset - cell array {dataset, context, celltype}
    %   mouseID - array indicating mouse ID for each dataset
    %   plot_info - structure with plotting parameters
    
    figure(700);clf
    positions = utils.calculateFigurePositions(1, 5, .5, []);
    if size( mod_index_by_dataset,2) > 2
        num_contexts = 2;
    else
        num_contexts = size( mod_index_by_dataset,2);
    end
    unique_mice = unique(mouseID);
    n_mice = length(unique_mice);
    n_celltypes =  size( mod_index_by_dataset,3);
    count = 0;
    x_lines = 0:num_contexts*n_celltypes+1;
    
    if nargin > 4
        ylims = varargin{1,1};
    end
    abs_logic = 1;
    if nargin > 5
        abs_logic = 0;
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
                    if abs_logic == 1
                        curr_data = abs(mod_index_by_dataset{d,context,celltype});
                    else
                        curr_data = (mod_index_by_dataset{d,context,celltype});
                    end
                    mouse_data = [mouse_data; curr_data(:)];
                end
                mouse_means(context) = nanmean(mouse_data);
            end
            
            % Plot connected line for this mouse
            if num_contexts > 1
                plot([x_pos(1)+.2,x_pos(2)-.2], mouse_means, '-', 'Color', ...
                    [plot_info.colors_celltypes(celltype,:), 0.3], 'LineWidth', 1)
            end
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
                mod_stats.stats(celltype,context).valid_datasets = find(~isnan(curr_means));
                 % Calculate 95% CI using bootstrapping
                mod_stats.stats(celltype,context).basic_stats =  get_basic_stats(valid_data);
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

            % One-sample permutation test vs zero
            n_perm = 10000;
            observed_mean = mean(valid_data);
            null_distribution = zeros(1, n_perm);
            for p = 1:n_perm
                signs = randi([0 1], size(valid_data)) * 2 - 1;  % random +1/-1
                null_distribution(p) = mean(valid_data .* signs);
            end
            p_val_vs_zero = mean(abs(null_distribution) >= abs(observed_mean));
            if p_val_vs_zero == 0
                p_val_vs_zero = 1/n_perm;
            end
            mod_stats.stats(celltype,context).p_val_vs_zero = p_val_vs_zero;

            if isfield(plot_info,'zero_star') && plot_info.zero_star == 1%only plot if requested
                if p_val_vs_zero < 0.05/(n_celltypes*num_contexts) %correct for multiple comparisons
%                     text(x_pos(context), mean_cel + err + 0.02, '*', ...
%                         'Color', plot_info.colors_celltypes(celltype,:), ...
%                         'FontSize', 8, 'HorizontalAlignment', 'center');
                            utils.plot_pval_star(x_pos(context),mean_cel + err + 0.02,p_val_vs_zero,[0,0],0,[0,0,0]);%plot_info.colors_celltypes(celltype,:)
                end
            end
        end

        % Statistical testing for this cell type
        if num_contexts > 1
            ct = 0;
            possible_tests = nchoosek(1:num_contexts,2);
            if max(mean_cell_all) > 0.1
                y_val = max(mean_cell_all) + 0.03;
            else
                y_val = max(mean_cell_all);
            end
        
            for t = 1:size(possible_tests,1)
                ctx1 = possible_tests(t,1);
                    ctx2 = possible_tests(t,2);
    
                    % Get data for statistical test
                    % Skip if this celltype has no data for either context
                    if isempty(mod_stats.stats) || ...
                       celltype > size(mod_stats.stats,1) || ...
                       isempty(mod_stats.stats(celltype,ctx1).valid_means) || ...
                       isempty(mod_stats.stats(celltype,ctx2).valid_means)
                        continue;
                    end
    
                    % Extract data
                    data1 = mod_stats.stats(celltype,ctx1).valid_means;
                    data2 = mod_stats.stats(celltype,ctx2).valid_means;
    
                    % Align datasets if lengths differ (use intersection of valid datasets if available)
                    if length(data1) ~= length(data2)
                        if isfield(mod_stats.stats(celltype,ctx1), 'valid_datasets') && ...
                           isfield(mod_stats.stats(celltype,ctx2), 'valid_datasets')
                           
                           % Intersect dataset indices from both contexts
                           common_datasets = intersect(mod_stats.stats(celltype,ctx1).valid_datasets, ...
                                                       mod_stats.stats(celltype,ctx2).valid_datasets);
                
                           [~, idx1] = ismember(common_datasets, mod_stats.stats(celltype,ctx1).valid_datasets);
                           [~, idx2] = ismember(common_datasets, mod_stats.stats(celltype,ctx2).valid_datasets);
                
                           data1 = data1(idx1);
                           data2 = data2(idx2);
                        else
                           % Skip if we can't align data
                           continue;
                        end
                    end
                
                if ~isempty(data1)
                [p_val_mod(t,celltype), ~, effectsize(t,celltype)] = permutationTest_updatedcb(...
                    data1, data2, 10000, 'paired', 1);
                else
                    p_val_mod(t,celltype) = 1;
                end
    
                if size(possible_tests,1) == 1 && size(mod_stats.stats,1) == 1 && celltype > 1%assume only pyr has valid stuff
                        p_val_mod(t,celltype) = 1;
                end
                
                if p_val_mod(t,celltype) < 0.05/n_celltypes
                    xline_vars = possible_tests(t,:) + ((celltype-1)*num_contexts);
                    ct = ct + 0.03;
                    utils.plot_pval_star(0, y_val+ct, p_val_mod(t,celltype), xline_vars, ...
                        0.01, plot_info.colors_celltypes(celltype,:))
                end
            end

        end
    end
    
    % Format plot
    xlim([x_lines(1) x_lines(end)])
    xticks(x_lines(2:end-1))
    xticklabels(repmat(plot_info.behavioral_contexts, 1, n_celltypes))

    if abs_logic == 1;
        ylabel({'Absolute Modulation';'Index'})
    else
        ylabel({'Modulation Index'})
    end

    if nargin > 6
        ylabel(varargin{1,3});
    end
    
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
%     set(gcf,'units','points','position',[10,100,(400/n_celltypes*num_contexts),170]);
    set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
    utils.set_current_fig;
    
    % Store statistics
    if num_contexts > 1
        mod_stats.tests = possible_tests;
        mod_stats.p_test = 'paired permutation across mice';
        mod_stats.p_val_mod = p_val_mod;
        mod_stats.effectsize = effectsize;
        mod_stats.n_mice = n_mice;
    end
    
    % Save results
    if ~isempty(save_dir)
        mkdir(save_dir)
        cd(save_dir)
        
        if nargin > 6
            save_string = varargin{1,3};
            save_string = strrep(save_string, '\', '');
            save_string = strrep(save_string, '/', '');
            save_string = strrep(save_string, '(', '');
            save_string = strrep(save_string, ')', '');
            saveas(700,strcat(save_string,'_abs',num2str(abs_logic), '_mod_index_connected_lines_n',num2str(n_mice),'.svg'));
            saveas(700,strcat(save_string,'_abs',num2str(abs_logic), '_mod_index_connected_lines_n',num2str(n_mice),'.fig'));
            exportgraphics(figure(700),strcat(save_string,'_abs',num2str(abs_logic), '_mod_index_connected_lines_n',num2str(n_mice),'_datasets.pdf'), 'ContentType', 'vector');
            save(strcat(save_string,'_abs',num2str(abs_logic), '_mod_index_stats_connected_lines_n',num2str(n_mice)),'mod_stats');

        else
            saveas(700,strcat('abs',num2str(abs_logic), '_mod_index_connected_lines_n',num2str(n_mice),'.svg'));
            saveas(700,strcat('abs',num2str(abs_logic), '_mod_index_connected_lines_n',num2str(n_mice),'.fig'));
            exportgraphics(figure(700),strcat('abs',num2str(abs_logic), '_mod_index_connected_lines_n',num2str(n_mice),'_datasets.pdf'), 'ContentType', 'vector');
            save(strcat('abs',num2str(abs_logic), '_mod_index_stats_connected_lines_n',num2str(n_mice)),'mod_stats');
        end
    end
end