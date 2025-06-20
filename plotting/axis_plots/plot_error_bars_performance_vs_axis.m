function binned_perf_all = plot_error_bars_performance_vs_axis(chosen_datasets,proj, axis_type, celltype,percent_correct,frame_range,edges_values,num_bins,save_dir)

% Bin engagement (trial_means) values
edges = linspace(edges_values(1), edges_values(2), num_bins+1); %linspace(min(all_proj), max(all_proj), num_bins+1);
bin_centers = edges(1:end-1) + diff(edges)/2;

% Plot
positions = utils.calculateFigurePositions(1, 5, .5, []);

ctx = 1; %focus on performance in the active context

num_datasets = length(chosen_datasets);
figure(802);clf;
hold on;


%
% performance across datasets
binned_perf_all = nan(num_bins, num_datasets); % num_bins x num_datasets

for dataset = chosen_datasets
    data = proj{dataset, celltype, ctx}.(axis_type);
    trial_means = mean(data(:, frame_range), 2);
    correct = percent_correct{dataset};

    for b = 1:num_bins
        bin_idx = trial_means >= edges(b) & trial_means < edges(b+1);
        trials_in_bin = correct(bin_idx);

        if ~isempty(trials_in_bin)
            binned_perf_all(b, dataset) = mean(trials_in_bin);
        end
    end
end

binned_perf = nanmean(binned_perf_all, 2);
binned_sem = nanstd(binned_perf_all, 0, 2) ./ sqrt(sum(~isnan(binned_perf_all), 2));

e = errorbar(bin_centers, binned_perf, binned_sem, '-o', 'LineWidth', 1,'MarkerSize',2,'Color','k');
set(e, 'CapSize', 2);   % Increase cap size (default is 6)
xlabel('Prestimulus "Engagement"');
ylabel({'Behavioral Performance';'(Fraction Correct)'});
ylim([0.5 1]);
xli = xlim;
xlim([xli(1) - .5,xli(2) + .5]); %adjust axis

set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

% Save results
if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(802,strcat('errorbar_correct_vs_axis_',num2str(axis_type),'_n_',num2str(length(chosen_datasets)),'.svg'));
    saveas(802,strcat('errorbar_correct_vs_axis_',num2str(axis_type),'_n_',num2str(length(chosen_datasets)),'.fig'));
    exportgraphics(figure(802),strcat('errorbar_correct_vs_axis_',num2str(axis_type),'_n_',num2str(length(chosen_datasets)),'.pdf'), 'ContentType', 'vector');
end