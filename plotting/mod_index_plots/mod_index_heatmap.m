function mod_index_heatmap(save_dir,stim_mod,plot_info,chosen_mice,varargin)
% Plots a heatmap of modulation indices across behavioral contexts.
%
% Inputs:
%   save_dir   - Directory to save figures (string, can be empty).
%   stim_mod   - [neurons x contexts] matrix of modulation indices.
%   plot_info  - Struct with field 'behavioral_contexts' (cell array of labels).
%   chosen_mice - Vector or cell array of mouse IDs used for labeling/saving.
%   varargin   - Optional:
%                   {1}: color axis limits, e.g. [-0.4 0.4]
%                   {2}: custom neuron sorting order (vector of indices)

figure(94);clf;
% cd 'C:\Code\Github\+colormaps'
% colorList= (colormaps.slanCM('coolwarm',100));
% colormap(colorList) % redblue
positions = utils.calculateFigurePositions(1, 5, .5, []);
mod_to_plot = stim_mod(:,1:length(plot_info.behavioral_contexts));
if 1 == length(plot_info.behavioral_contexts)
    %halve positions
    positions(1,3) = positions(1,3) * 0.2;  % width is 3rd element
end
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
if length(plot_info.behavioral_contexts) > 1
xticklabels(plot_info.behavioral_contexts);
else
    xticklabels('');
end
ylabel('Neurons','FontSize',7)
set(heatmap_plot,'AlphaData',~isnan(mod_to_plot(id,:))); %,~isnan(stim_mod(chosen_cels(id),:)))
set(gca,'color',[.7 .7 .7],'fontsize',7);
% set(gcf,'position',[100,100,200,200])
set(gca, 'FontSize', 7, 'Units', 'inches', 'Position', positions(1, :));
movegui(gcf, 'center')

if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(figure(94),strcat('mod_index_heatmap_',num2str(length(chosen_mice)),'_datasets.svg'));
    saveas(figure(94),strcat('mod_index_heatmap_',num2str(length(chosen_mice)),'_datasets.fig'));
    exportgraphics(figure(94),strcat('mod_index_heatmap_',num2str(length(chosen_mice)),'_datasets.pdf'), 'ContentType', 'vector');
    save('sorting_id_heatmap','id');
end