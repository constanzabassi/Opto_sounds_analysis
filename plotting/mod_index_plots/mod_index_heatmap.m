function mod_index_heatmap(save_dir,stim_mod,plot_info,chosen_mice,varargin)
figure(94);clf;
% cd 'C:\Code\Github\+colormaps'
% colorList= (colormaps.slanCM('coolwarm',100));
% colormap(colorList) % redblue
mod_to_plot = stim_mod(:,1:length(plot_info.behavioral_contexts));
colormap redblue
meancontextmod = nanmean(mod_to_plot(:,:) , 2);%nanmean(stim_mod(chosen_cels,:) , 2);

[~ ,id] = sort(meancontextmod,'descend');
if nargin > 5
     id = varargin{1,2};
end
heatmap_plot =imagesc(mod_to_plot(id,:)); %imagesc(stim_mod(chosen_cels(id),:));

%plot(run_mod')
if nargin > 4
     caxis(varargin{1,1})
else
    caxis([ -.4 .4])
end
colorbar
xticks([1:size(mod_to_plot,2)]);
xticklabels(plot_info.behavioral_contexts);
ylabel('Neurons')
set(heatmap_plot,'AlphaData',~isnan(mod_to_plot(id,:))); %,~isnan(stim_mod(chosen_cels(id),:)))
set(gca,'color',[.7 .7 .7],'fontsize',12);
set(gcf,'position',[100,100,200,200])
movegui(gcf, 'center')

% axis square
% box off
% set(gcf,'position',[100,100,300,180])
% movegui(gcf, 'center')
% utils.set_current_fig;

if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(figure(94),strcat('mod_index_heatmap_',num2str(length(chosen_mice)),'_datasets.svg'));
    saveas(figure(94),strcat('mod_index_heatmap_',num2str(length(chosen_mice)),'_datasets.fig'));
    exportgraphics(figure(94),strcat('mod_index_heatmap_',num2str(length(chosen_mice)),'_datasets.pdf'), 'ContentType', 'vector');
    save('sorting_id_heatmap','id');
end