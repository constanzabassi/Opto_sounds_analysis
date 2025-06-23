function [binned_resp_all_celltype,errorbar_resp_ct_stats] = plot_error_bars_response_celltypes_vs_axis(chosen_datasets,proj, axis_type,proj2, axis_type2, ctx,frame_range1,frame_range2,edge_values,num_bins,colorss,save_dir)
num_datasets = length(chosen_datasets);
positions = utils.calculateFigurePositions(1, 5, .5, []);

% Bin by engagement
edges = linspace(edge_values(1),edge_values(2), num_bins+1); %linspace(min(all_engagement), max(all_engagement), num_bins+1);
bin_centers = edges(1:end-1) + diff(edges)/2;

if ctx == 1
    string_context = 'Active';
else
    string_context = 'Passive';
end
figure(803);clf;
hold on;
for celltype = 1:3
    % performance across datasets
    binned_resp_all = nan(num_bins, num_datasets); % num_bins x num_datasets
    
    for dataset = chosen_datasets
        if ~isempty(proj{dataset, celltype, ctx})
            data = proj{dataset, celltype, ctx}.(lower(axis_type));
            trial_means = mean(data(:, frame_range1), 2);
        else
            trial_means =nan;
        end
        
        if ~isempty(proj2{dataset, celltype, ctx})
            data2 = proj2{dataset, celltype, ctx}.(lower(axis_type2));
            trial_means2 = mean(data2(:, frame_range2), 2);
        else
            trial_means2 =nan;
        end
    
        for b = 1:num_bins
            bin_idx = trial_means >= edges(b) & trial_means < edges(b+1);
            values_in_bin = trial_means2(bin_idx);
    
            if ~isempty(values_in_bin )
                binned_resp_all(b, dataset) = mean(values_in_bin);
            end
        end
    end
    binned_resp_all_celltype{celltype} = binned_resp_all;
    binned_resp{celltype} = nanmean(binned_resp_all, 2);
    binned_sem = nanstd(binned_resp_all, 0, 2) ./ sqrt(sum(~isnan(binned_resp_all), 2));
    
    % Plot
    e = errorbar(bin_centers, binned_resp{celltype}, binned_sem, '-o','LineWidth', 1,'MarkerSize',2,'Color',colorss(celltype,:));
    set(e, 'CapSize', 2);   % Increase cap size (default is 6)

end
xlabel('Prestimulus "Engagement"');
ylabel({[string_context, ' ', axis_type2]; 'Neural Response'});
xli = xlim;
xlim([xli(1) - .5,xli(2) + .5]); %adjust axis

%Do statistical comparisons across bins
celltype_combos = nchoosek(1:3, 2); % All 3 pairwise combinations
num_combos = size(celltype_combos, 1);
p_values = nan(num_bins, num_combos);

for b = 1:num_bins
    for c = 1:num_combos
        ct1 = celltype_combos(c, 1);
        ct2 = celltype_combos(c, 2);

        x = binned_resp_all_celltype{ct1}(b, :);
        y = binned_resp_all_celltype{ct2}(b, :);
        valid = ~isnan(x) & ~isnan(y);

        errorbar_resp_ct_stats.stats{b, ct1} = utils.get_basic_stats(x);
        errorbar_resp_ct_stats.stats{b, ct2} = utils.get_basic_stats(y);

        if any(valid)
            p_values(b, c) = permutationTest(x(valid), y(valid), 10000); % unpaired
        end
    end
end

% Store results
errorbar_resp_ct_stats.p_values = p_values;
errorbar_resp_ct_stats.test = 'unpaired permutation across bins and celltype pairs';
errorbar_resp_ct_stats.combos = celltype_combos;
errorbar_resp_ct_stats.significant = find(p_values < 0.05 / (num_bins * num_combos)); % Bonferroni correction

%add significant stars
if isfield(errorbar_resp_ct_stats, 'p_values') && ~isempty(errorbar_resp_ct_stats.p_values)
    num_bins = size(errorbar_resp_ct_stats.p_values, 1);
    ct_combos = nchoosek(1:3, 2); % assuming 3 cell types
    offset = 0.05; % vertical spacing between stars

    for b = 1:num_bins
        for c = 1:size(ct_combos, 1)
            p = errorbar_resp_ct_stats.p_values(b, c);
            if p < 0.05 / num_bins  % Bonferroni correction
                y = ylim;
                star_y = y(2) - 0.05 + (c * offset);
                x_pos = bin_centers(b);
                utils.plot_pval_star(0, star_y, p, [x_pos, x_pos], 0.05);
            end
        end
    end
end



set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

% Save results
if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
%     saveas(803,strcat('errorbar_response_celltypes_vs_axis_',num2str(axis_type2),'_ctx_',num2str(ctx),'_n_',num2str(length(chosen_datasets)),'_edges_',num2str(edge_values),'.svg'));
    saveas(803,strcat('errorbar_response_celltypes_vs_axis_',num2str(axis_type2),'_ctx_',num2str(ctx),'_n_',num2str(length(chosen_datasets)),'_edges_',num2str(edge_values),'.fig'));
    exportgraphics(figure(803),strcat('errorbar_response_celltypes_vs_axis_',num2str(axis_type2),'_ctx_',num2str(ctx),'_n_',num2str(length(chosen_datasets)),'_edges_',num2str(edge_values),'.pdf'), 'ContentType', 'vector');

    save(strcat('errorbar_resp_celltypes_vs_axis_',num2str(axis_type2),'_ctx_',num2str(ctx),'_stats_n',num2str(length(chosen_datasets)),'_edges_',num2str(edge_values)),'errorbar_resp_ct_stats');

end
