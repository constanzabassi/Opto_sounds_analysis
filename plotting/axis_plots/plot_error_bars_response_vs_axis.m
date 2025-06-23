function [binned_resp_all_ctx,errorbar_resp_stats] = plot_error_bars_response_vs_axis(chosen_datasets,proj, axis_type,proj2, axis_type2, celltype,frame_range1,frame_range2,edge_values,num_bins,colorss,save_dir)
num_datasets = length(chosen_datasets);
positions = utils.calculateFigurePositions(1, 5, .5, []);

% Bin by engagement
edges = linspace(edge_values(1),edge_values(2), num_bins+1); %linspace(min(all_engagement), max(all_engagement), num_bins+1);
bin_centers = edges(1:end-1) + diff(edges)/2;


figure(803);clf;
hold on;
for ctx = 1:2
    % performance across datasets
    binned_resp_all = nan(num_bins, num_datasets); % num_bins x num_datasets
    
    for dataset = chosen_datasets
        data = proj{dataset, celltype, ctx}.(lower(axis_type));
        trial_means = mean(data(:, frame_range1), 2);
        
        data2 = proj2{dataset, celltype, ctx}.(lower(axis_type2));
        trial_means2 = mean(data2(:, frame_range2), 2);
    
        for b = 1:num_bins
            bin_idx = trial_means >= edges(b) & trial_means < edges(b+1);
            values_in_bin = trial_means2(bin_idx);
    
            if ~isempty(values_in_bin )
                binned_resp_all(b, dataset) = mean(values_in_bin);
            end
        end
    end
    binned_resp_all_ctx{ctx} = binned_resp_all;
    binned_resp{ctx} = nanmean(binned_resp_all, 2);
    binned_sem = nanstd(binned_resp_all, 0, 2) ./ sqrt(sum(~isnan(binned_resp_all), 2));
    
    % Plot
    e = errorbar(bin_centers, binned_resp{ctx}, binned_sem, '-o','LineWidth', 1,'MarkerSize',2,'Color',colorss(ctx,:));
    set(e, 'CapSize', 2);   % Increase cap size (default is 6)

end
xlabel('Prestimulus "Engagement"');
ylabel({strcat(axis_type2, ' Neural'); 'Response'});
xli = xlim;
xlim([xli(1) - .5,xli(2) + .5]); %adjust axis

%Do statistical comparisons across bins
p_values = nan(num_bins, 1);
for b = 1:(num_bins)
    x = binned_resp_all_ctx{1}(b, :);
    y = binned_resp_all_ctx{2}(b, :);
    valid = ~isnan(x) & ~isnan(y);
    errorbar_resp_stats.stats{b,1} = utils.get_basic_stats(x);
    errorbar_resp_stats.stats{b,2} = utils.get_basic_stats(y);

    errorbar_resp_stats.p_values(b) = permutationTest_updatedcb(x(valid), y(valid), 10000, 'paired', 1);
end
errorbar_resp_stats.test = 'paired permutation across bins';
errorbar_resp_stats.significant = find(errorbar_resp_stats.p_values < 0.05/num_bins);

%add significant stars
ct = 0;
if ~isempty(errorbar_resp_stats.significant)
    for i = errorbar_resp_stats.significant
        yli = ylim;
        xli = xlim;
        utils.plot_pval_star(0, yli(2)-(.025*length(errorbar_resp_stats.significant))+ct, errorbar_resp_stats.p_values(i),[bin_centers(i),bin_centers(i)],.05);
        ct = ct+.07;
    end
end


set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

legend_string = {'Active','Passive'};
% Get current axis limits
x_range = xlim;
y_range = ylim;
% Calculate base text position
text_x = x_range(2) -.09 * diff(x_range);
text_y = y_range(2) - .2 * diff(y_range);

% Auto-calculate evenly spaced y-offsets
num_labels = 2;
y_offsets = linspace(0, 0.1 * (num_labels - 1), num_labels); % Adjusted scaling
% Place text labels
for i = 1:num_labels
    text(text_x, text_y - y_offsets(i) * diff(y_range), legend_string{i}, ...
         'Color',colorss(i,:), 'FontSize', 8);
end


% Save results
if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(803,strcat('errorbar_response_vs_axis_',num2str(axis_type2),'_n_',num2str(length(chosen_datasets)),'.svg'));
    saveas(803,strcat('errorbar_response_vs_axis_',num2str(axis_type2),'_n_',num2str(length(chosen_datasets)),'.fig'));
    exportgraphics(figure(803),strcat('errorbar_response_vs_axis_',num2str(axis_type2),'_n_',num2str(length(chosen_datasets)),'.pdf'), 'ContentType', 'vector');

    save(strcat('errorbar_resp_vs_axis_',num2str(axis_type2),'_stats_n',num2str(length(chosen_datasets))),'errorbar_resp_stats');

end
