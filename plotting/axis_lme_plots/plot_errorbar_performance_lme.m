function lme_perf = plot_errorbar_performance_lme(percent_correct_concat,engagement_proj_all_sound,engagement_proj_all_stim,context_all_sound,context_all_stim,save_dir)
%% collect correctness across mice
correct_all_ctrl = []; animal_id_all=[];
for dataset = 1:24
    correct_all_ctrl = [correct_all_ctrl; percent_correct_concat{dataset}'];
    animal_id_all = [animal_id_all; repmat(dataset, size(percent_correct_concat{dataset},2), 1)];
end

active_trials = context_all_sound == 0;
engagement_active = [engagement_proj_all_sound(active_trials);engagement_proj_all_stim(context_all_stim == 0)];
performance_active =  correct_all_ctrl;
animal_active = animal_id_all;  % Categorical ID

%concatenate with stim trials to have more trials!

%% me logistic regression

% Create table
tbl_perf = table(engagement_active, performance_active, animal_active, ...
    'VariableNames', {'EngagementProj', 'Performance', 'AnimalID'});
tbl_perf.Performance = logical(tbl_perf.Performance); % Ensure binary
tbl_perf.AnimalID = categorical(tbl_perf.AnimalID);   % Random effect grouping

% Fit GLME: random intercept for each animal
lme_perf = fitglme(tbl_perf, 'Performance ~ EngagementProj + (1|AnimalID)', ...
                   'Distribution', 'binomial');
disp(lme_perf)


edges = prctile(engagement_active, linspace(0, 100, 6));  % 5 bins
bin_centers = movmean(edges, 2, 'Endpoints','discard');
n_bins = length(bin_centers);

% Initialize
% animals = categories(animal_id_all);
n_animals = 24;
perf_by_bin = nan(n_animals, n_bins);

% Loop through animals
for i = 1:n_animals
    idx = animal_id_all == i;
    
    e = engagement_active(idx);
    p = performance_active(idx);

    if numel(e) > 10
        [~,~,bin] = histcounts(e, edges);
        for b = 1:n_bins
            trials_in_bin = p(bin == b);
            if numel(trials_in_bin) >= 3
                perf_by_bin(i, b) = mean(trials_in_bin);
            end
        end
    end
end

% Compute mean and SEM across animals
mean_perf = nanmean(perf_by_bin, 1);
sem_perf = nanstd(perf_by_bin, 0, 1) ./ sqrt(sum(~isnan(perf_by_bin), 1));

% Plot
figure(105);clf; hold on;
errorbar(bin_centers, mean_perf, sem_perf, '-ok', 'MarkerFaceColor', 'k', 'LineWidth', 1,'MarkerSize',3,'CapSize',2);
xlabel('Engagement projection');
ylabel('Fraction Correct');
% title('Performance vs Engagement (mean ± SEM across mice)');
set(gca, 'FontSize', 8);
positions = utils.calculateFigurePositions(1, 5, .5, []);
set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;
xli = xlim;
xlim([xli(1)- xli(2)*.3,xli(2) + xli(2)*.3]); %adjust axis
ylim([.6 .85])

figure(106);clf
bar(sum(~isnan(perf_by_bin), 1));
xlabel('Bin'); ylabel('# Animals'); title('Animals per Engagement Bin');

if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(105,strcat('errorbar_performance_vs_engagement_axis','.fig'));
    exportgraphics(figure(105),strcat('errorbar_performance_vs_engagement_axis','.pdf'), 'ContentType', 'vector');
end
