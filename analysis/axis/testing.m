function plot_projected_trials
figure;
hold on;
for contexts = 1:2
    plot(1:122,mean(sound_proj_ctrl(total_trials{current_dataset, contexts, 2},:)),'color', colorss(contexts,:)); % 
end

%%
current_dataset = 1;
celltype = 4;
colorss = [0,0,0;.5,.5,.5];
figure;
hold on;
for dataset = 1:24
for contexts = 1:2
    plot(1:122,(mean(proj_ctrl{dataset,celltype,contexts}.sound) - mean(proj_ctrl{dataset,celltype,contexts}.sound(:,1:59),'all')),'color', colorss(contexts,:)); % 
end
end

figure;
hold on;
for dataset = 1:24
for contexts = 1:2
    plot(1:122,(mean(proj_ctrl{dataset,celltype,contexts}.sound)),'color', colorss(contexts,:)); % 
end
end

figure;
hold on;
for dataset = 1:24
for contexts = 1:2
    plot(1:122,(mean(proj{dataset,celltype,contexts}.stim) - mean(proj{dataset,celltype,contexts}.stim(:,1:59),'all')),'color', colorss(contexts,:)); % 
end
end

% figure;
% hold on;
% for dataset = 1:24
% for contexts = 1:2
%     plot(1:122,(mean(proj_ctrl{dataset,celltype,contexts}.context) - mean(proj_ctrl{dataset,celltype,contexts}.context(:,61:122),'all')),'color', colorss(contexts,:)); % 
% end
% end

figure;
hold on;
for dataset = 1:24
for contexts = 1:2
    plot(1:122,proj_ctrl{dataset,celltype,contexts}.context,'color', colorss(contexts,:)); % 
end
end


temp = [];temp2 =[];
figure;
hold on;
for dataset = 1:size(trial_corr_context,1)
    scatter(trial_corr_context{dataset,celltype,1}.stim,trial_corr_context{dataset,celltype,2}.stim,'magenta'); % 
    scatter(trial_corr_context{dataset,celltype,1}.sound,trial_corr_context{dataset,celltype,2}.sound,'cyan'); % 
    temp = [temp;trial_corr_context{dataset,celltype,1}.sound,trial_corr_context{dataset,celltype,2}.sound];
    temp2 = [temp2;trial_corr_context{dataset,celltype,1}.stim,trial_corr_context{dataset,celltype,2}.stim];
end
%%
figure;
hold on;

% Initialize
num_timepoints = 122;
num_datasets = 24;
contexts = 2;

% Preallocate
all_data = cell(contexts, 1);

for ctx = 1:contexts
    ctx_data = zeros(num_datasets, num_timepoints);
    for dataset = 1:num_datasets
        data = proj_ctrl{dataset, celltype, ctx}.sound; %proj{dataset,celltype,ctx}.stim;%proj_ctrl{dataset, celltype, ctx}.sound; %proj{dataset,celltype,ctx}.stim;%proj{dataset,celltype,ctx}.stim%proj_ctrl{dataset, celltype, ctx}.sound;
        baseline = 0; %mean(data(:,1:59), 'all');
        ctx_data(dataset, :) = mean(data) - baseline;
    end
    all_data{ctx} = ctx_data;
end

% Plot with mean and SEM
for ctx = 1:contexts
    data = all_data{ctx};
    m = mean(data, 1);
    s = std(data, [], 1) / sqrt(size(data,1)); % SEM

    % Plot shaded error bar or mean with error bars
    fill([1:num_timepoints, fliplr(1:num_timepoints)], ...
         [m + s, fliplr(m - s)], ...
         colorss(ctx,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot(1:num_timepoints, m, 'Color', colorss(ctx,:), 'LineWidth', 2);
end

xlabel('Time');
ylabel('Mean-centered projection');
legend({'Context 1', 'Context 2'});
title('Mean ± SEM of Sound Projection by Context');

%%
figure;
hold on;

frame_range = 50:60; %63:93; % specify the frames you're averaging over
bin_edges = -1.5:0.2:1.5; % adjust bins as needed
all_trial_means_contexts = {};
all_trial_means_datasets = [];
for ctx = 1:2
    all_trial_means = [];

    for dataset = 1:24
        data = proj_ctrl{dataset, celltype, ctx}.context; %proj{dataset,celltype,ctx}.stim; %proj_ctrl{dataset, celltype, ctx}.sound;
        trial_means = mean(data(:, frame_range), 2); % mean per trial in specified frames
        all_trial_means = [all_trial_means; trial_means];
        all_trial_means_datasets(ctx,dataset) = mean(trial_means);
    end

    % Compute histogram (normalized to get fraction of trials)
    h = histogram(all_trial_means, bin_edges, 'Normalization', 'probability', ...
                  'DisplayStyle', 'stairs', 'LineWidth', 2, 'EdgeColor', colorss(ctx,:));
    all_trial_means_contexts{ctx} = all_trial_means;

end

xlabel('Mean projection (frames 60–90)');
ylabel('Fraction of trials');
legend({'Context 1', 'Context 2'});
title('Histogram of Trial Means by Context');

num_contexts = 2;
possible_tests = nchoosek(1:num_contexts,2);
for t = 1:size(possible_tests,1)
            ctx1 = possible_tests(t,1);
                ctx2 = possible_tests(t,2);
                data1 = all_trial_means_contexts{ctx1};
                data2 = all_trial_means_contexts{ctx2};
[p_val_mod(t,celltype), ~, effectsize(t,celltype)] = permutationTest(...
data1, data2, 10000);
end

num_contexts = 2;
possible_tests = nchoosek(1:num_contexts,2);
for t = 1:size(possible_tests,1)
            ctx1 = possible_tests(t,1);
                ctx2 = possible_tests(t,2);
                data1 = all_trial_means_datasets(ctx1,:);
                data2 = all_trial_means_datasets(ctx2,:);
[p_val_mod2(t,celltype), ~, effectsize2(t,celltype)] = permutationTest_updatedcb(...
data1, data2, 10000,'paired',1);
end


