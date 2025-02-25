function make_heatmap_sorted(data,x_label,y_label,ylims,varargin)

colormap redblue
data_to_plot= squeeze(data);

if nargin > 4
    imagesc(data_to_plot(varargin{1,1},:)); 
else
    imagesc(data_to_plot); 
end

clim(ylims);
% xlim([0 size(data,2)]);
% ylim([0 size(data,1)]);
if length(x_label) > 1
    xticks([1:size(data,2)]);
    xticklabels(x_label);
else
    xlabel(x_label);
end
ylabel(y_label);