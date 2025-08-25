function plot_avg_trial_traces_simple(data,sig_neurons,mod_index,dataset_index,plot_mode,save_dir)
if length(sig_neurons) > 5
        total_to_plot = 5;
    else
        total_to_plot = length(sig_neurons);
end
% Validate plot_mode.
valid_modes = {'both', 'stim', 'ctrl'};
if ~ismember(plot_mode, valid_modes)
    error('Unknown plot_mode. Use ''both'', ''stim'', or ''ctrl''.');
end

for i = 1:total_to_plot
    neuron_id = sig_neurons(i);
    figure(200); clf; 
    hold on
    for tr=1:size(data,1) %1:length trials
        plot(squeeze(data(tr,neuron_id,:)),'color',[0.7 0.7 0.7],'LineWidth',.3);
    end
    
    plot(squeeze(mean(data(:,neuron_id,:))),'k','LineWidth',1.0);
    % Hide axis ticks and labels, but keep title
    set(gca,'XTick',[],'YTick',[]);
    set(gca,'XColor','none','YColor','none');

    title(sprintf('Mod Index: %.2f', mod_index(neuron_id)), 'FontWeight', 'normal','Fontsize',6);
    
    set(gca, 'box', 'off')
    % set stim onset line
    if strcmp(plot_mode,'ctrl')
        xline(61, 'y','color',[0.5 0.5 0.5],'LineWidth',2,'alpha',0.5);% [0.9 0.8 0.2]
    else
        xline(61, 'y','color',[0.9290 0.6940 0.1250],'LineWidth',2,'alpha',0.5);% [0.9 0.8 0.2]
    end

    %DRAW SCALE BARS FOR X AND Y AXIS
    % Get current axis limits
    xlims = xlim;
    ylims = ylim;
    % Determine Y scale based on ylim
    if ylims(2) > 1
        y_scale = 1;
    else
        y_scale = 0.2;
    end
    % Scale bar lengths
    x_scale = 30; % 1 second
    % Position: bottom-left outside plot (slightly offset)
    offset_x = xlims(1) - 0.03 * range(xlims); % small shift left
    offset_y = ylims(1) - 0.03 * range(ylims); % small shift down
    hold on;
    % Draw X scale bar
    plot([offset_x offset_x + x_scale], [offset_y offset_y], 'k', 'LineWidth', 1);
    % Draw Y scale bar
    plot([offset_x offset_x], [offset_y offset_y + y_scale], 'k', 'LineWidth', 1);
    axis square
    set(gcf,'Position',[100,100,50,70])
    hold off
    movegui(gcf,'center');
    
    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end
        saveas(gcf, fullfile(save_dir, sprintf('Avg_traces_neuron_%d_dataset_%i_context_%c.fig', neuron_id, dataset_index(1), num2str(dataset_index(2)))));
        exportgraphics(gcf, fullfile(save_dir, sprintf('Avg_traces_neuron_%d_dataset_%i_context_%c.pdf', neuron_id, dataset_index(1), num2str(dataset_index(2)))), 'ContentType', 'vector');
    end
end
