function plot_corr_matrix_projections(proj,celltype,frames_pre,frames_post,save_dir)

all_sound = [];
all_eng = [];
all_stim = [];
for dataset = 1:size(proj,1)
    all_sound = [all_sound,mean(proj{dataset,celltype,1}.sound(:,frames_post),2)',mean(proj{dataset,celltype,2}.sound(:,frames_post),2)'];
    all_stim = [all_stim,mean(proj{dataset,celltype,1}.stim(:,frames_post),2)',mean(proj{dataset,celltype,2}.stim(:,frames_post),2)'];
    all_eng = [all_eng,mean(proj{dataset,celltype,1}.context(:,frames_pre),2)',mean(proj{dataset,celltype,2}.context(:,frames_pre),2)'];

end

proj_matrix = [all_eng; all_sound; all_stim]';  % [trials x 3]
[corr_mat,p_vals] = corrcoef(proj_matrix);

figure(801);clf;
colormap parula
imagesc(corr_mat);
colorbar;
xticks(1:3); xticklabels({'Engagement','Sound','Stim'});
yticks(1:3); yticklabels({'Engagement','Sound','Stim'});

positions = utils.calculateFigurePositions(1, 5, .5, []);
set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(801,strcat('proj_correlation_matrix.fig'));
    exportgraphics(figure(801),strcat('proj_correlation_matrix.pdf'), 'ContentType', 'vector');
%     save(strcat('hist_stats_n',num2str(length(chosen_datasets))),'hist_stats');
end