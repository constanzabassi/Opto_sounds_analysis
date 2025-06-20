function plot_error_bars_response_vs_axis(chosen_datasets,proj, axis_type,proj2, axis_type2, celltype,frame_range1,frame_range2,edge_values,num_bins,colorss,save_dir)
num_datasets = length(chosen_datasets);
positions = utils.calculateFigurePositions(1, 5, .5, []);

% Bin by engagement
edges = linspace(edge_values(1),edge_values(2), num_bins+1); %linspace(min(all_engagement), max(all_engagement), num_bins+1);
bin_centers = edges(1:end-1) + diff(edges)/2;


figure(802);clf;
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
    
    binned_perf = nanmean(binned_resp_all, 2);
    binned_sem = nanstd(binned_resp_all, 0, 2) ./ sqrt(sum(~isnan(binned_resp_all), 2));
    
    % Plot
    errorbar(bin_centers, binned_perf, binned_sem, '-o','LineWidth', 1,'MarkerSize',2,'Color',colorss(ctx,:));
end
xlabel('Prestimulus "Engagement"');
ylabel({strcat(axis_type2, ' Neural'); 'Response'});

set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

% Save results
if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(803,strcat('errorbar_response_vs_axis_',num2str(axis_type),'_n_',num2str(length(chosen_datasets)),'.svg'));
    saveas(803,strcat('errorbar_response_vs_axis_',num2str(axis_type),'_n_',num2str(length(chosen_datasets)),'.fig'));
    exportgraphics(figure(803),strcat('errorbar_correct_vs_axis_',num2str(axis_type),'_n_',num2str(length(chosen_datasets)),'.pdf'), 'ContentType', 'vector');
end
