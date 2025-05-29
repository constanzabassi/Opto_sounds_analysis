function plot_avg_traces_baseline_subtracted(deconv_response, colors, lineStyles_contexts, celltypes_ids, frames, stim_frame, save_dir, average_over_neurons,type,plot_info)

if nargin < 9
    average_over_neurons = false;
end

positions = utils.calculateFigurePositions(1, 6, .4, []);
contexts = {'active', 'passive'};
data_modes = {'raw', 'bs'}; % raw and baseline subtracted
stim_ctrl_idx = [1, 0, 1, 0];

for fig_idx = 1:4
    figure(fig_idx); clf;

    for celtype = 1:size(deconv_response,3)
        subplot(1,size(deconv_response, 3),celtype);
        hold on;

        for context = 1:size(deconv_response,1)
            all_traces = [];  % All neuron traces pooled
            mouse_means = []; % Mean trace per dataset

            for mouse = 1:size(deconv_response,2)
                dat_struct = deconv_response{context,mouse,celtype};
                if isempty(dat_struct)
                    continue
                end

                if stim_ctrl_idx(fig_idx) == 1
                    if ~isfield(dat_struct, 'stim') || isempty(dat_struct.stim)
                        continue
                    end
                    if length(size(dat_struct.stim))>2
                        dat = squeeze(mean(dat_struct.stim));
%                         if contains(type,'deconv')
%                             frame_window = 1:size(dat_struct.stim,3);
%                             dat = squeeze(squeeze(sum(dat_struct.stim(:,:,frame_window),1))/(length(frame_window)/30));
%                         end
                    else
                        dat = dat_struct.stim;
                    end
                else
                    if ~isfield(dat_struct, 'ctrl') || isempty(dat_struct.ctrl)
                        continue
                    end
                    if length(size(dat_struct.stim))>2
                        dat = squeeze(mean(dat_struct.ctrl));
                    else
                        dat = dat_struct.ctrl;
                    end

                end

                if isnan(dat)
                    continue
                end

                dat = dat(:, frames);

                if strcmp(data_modes{ceil(fig_idx/2)}, 'bs')
                    baseline = mean(dat(:, 31:60), 2);
                    dat = dat - baseline;
                end

                if average_over_neurons
                    all_traces = [all_traces; dat];  % Pool all neurons
                else
                    mean_trace = mean(dat, 1);  % Average per dataset
                    mouse_means = [mouse_means; mean_trace];
                end
            end

            if average_over_neurons && ~isempty(all_traces)
                shadedErrorBar([], mean(all_traces,1), ...
                    std(all_traces,[],1)/sqrt(size(all_traces,1)), ...
                    'lineprops', {'LineWidth', 1., 'LineStyle', '-', ...
                    'Color', colors((celtype-1)*3+context,:)});
            elseif ~isempty(mouse_means)
                shadedErrorBar([], mean(mouse_means,1), ...
                    std(mouse_means,[],1)/sqrt(size(mouse_means,1)), ...
                    'lineprops', {'LineWidth', 1., 'LineStyle', '-', ...
                    'Color', colors((celtype-1)*3+context,:)});
            end

%             %place legends (plot too small to do this consistently)
%             if context == size(deconv_response,1)
%                 utils.place_text_labels(plot_info.behavioral_contexts, colors((celtype-1)*3+1:(celtype-1)*3+2,:), 0.75, 8); %.4
%             end
                hold off;
        end

        xlimss = [31 91];
        xlim(xlimss );
        xticks([31 61 91]);
        xticklabels([-1 0 1]);

%         [xticks_in, xticks_lab] = utils.x_axis_sec_aligned(30, length(xlimss ));
%         xticks(xticks_in);
%         xticklabels(xticks_lab);

        title(celltypes_ids{celtype}, 'FontSize', 8, 'FontName', 'arial','FontWeight','normal');


        
        set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(celtype, :));
        
        yli = ylim;
        for f = 1:size(stim_frame,1)
            if contains(type,'opto')
                color_onset = [1 0.8 0.3];
            else
                color_onset = [.5 .5 .5];
            end

            rectangle('Position', [stim_frame(f,1), yli(1), stim_frame(f,2)-stim_frame(f,1), ...
                yli(2)-yli(1)], 'FaceColor', color_onset, 'EdgeColor', 'none');
        end
%         xline(stim_frame(1), '--k', 'LineWidth', 1);

    end

    if ~isempty(save_dir)
        fig_suffix = {'raw_stim','raw_ctrl','bs_stim','bs_ctrl'};
        mkdir(fullfile(save_dir, 'avg_traces'));
        saveas(fig_idx, fullfile(save_dir, 'avg_traces', ...
            strcat('avg_traces_', fig_suffix{fig_idx},'_', type,'.fig')));
        exportgraphics(gcf, fullfile(save_dir, 'avg_traces', ...
            strcat('avg_traces_', fig_suffix{fig_idx},'_', type,'.pdf')), 'ContentType', 'vector');
    end
end

% 
% positions = utils.calculateFigurePositions(1, 5, .5, []);
% contexts = {'active', 'passive'};
% data_modes = {'raw', 'bs'}; % raw and baseline subtracted
% stim_ctrl_idx = [1, 0, 1, 0];
% 
% for fig_idx = 1:4
%     figure(fig_idx); clf;
% %     t = tiledlayout(1, size(deconv_response, 3), 'TileSpacing', 'Compact', 'Padding', 'Compact');
% 
%     for celtype = 1:size(deconv_response,3)
% %         nexttile;
%         subplot(1,size(deconv_response, 3),celtype);
%         hold on;
% 
%         for context = 1:2
%             mouse_means = [];
% 
%             for mouse = 1:size(deconv_response,2)
%                 dat_struct = deconv_response{context,mouse,celtype};
%                 if isempty(dat_struct)
%                     continue
%                 end
% 
%                 if stim_ctrl_idx(fig_idx) == 1
%                     if ~isfield(dat_struct, 'stim') || isempty(dat_struct.stim)
%                         continue
%                     end
%                     dat = squeeze(mean(dat_struct.stim));
%                 else
%                     if ~isfield(dat_struct, 'ctrl') || isempty(dat_struct.ctrl)
%                         continue
%                     end
%                     dat = squeeze(mean(dat_struct.ctrl));
%                 end
% 
%                 if isnan(dat)
%                     continue
%                 end
% 
%                 dat = dat(:, frames);
% 
%                 if strcmp(data_modes{ceil(fig_idx/2)}, 'bs')
%                     baseline = mean(dat(:, 31:60), 2);
%                     dat = dat - baseline;
%                 end
% 
%                 mean_trace = mean(dat, 1);
%                 mouse_means = [mouse_means; mean_trace];
%             end
% 
%             if ~isempty(mouse_means)
%                 shadedErrorBar([], mean(mouse_means,1), ...
%                     std(mouse_means,[],1)/sqrt(size(mouse_means,1)), ...
%                     'lineprops', {'LineWidth', 1.5, 'LineStyle', lineStyles_contexts{context}, 'Color', colors((celtype-1)*3+context,:)});
%             end
%         end
% 
%         xlim([31 91]);
%         xticks([31 61 91]);
%         xticklabels([-1 0 1]);
%         [xticks_in, xticks_lab] = utils.x_axis_sec_aligned(stim_frame(1,1), length([31:91]));
%         xticks(xticks_in);
%         xticklabels(xticks_lab);
% 
%         title(celltypes_ids{celtype}, 'FontSize', 8, 'FontName', 'arial','FontWeight','normal');
% 
%         yli = ylim;
%         for f = 1:size(stim_frame,1)
%             rectangle('Position', [stim_frame(f,1), yli(1), stim_frame(f,2)-stim_frame(f,1), yli(2)-yli(1)], ...
%                       'FaceColor', [1 0.8 0.3], 'EdgeColor', 'none');
%         end
%         xline(stim_frame(1), '--k', 'LineWidth', 1);
%         hold off;
%         set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(celtype, :));
%     end
%     
% 
%     if ~isempty(save_dir)
%         fig_suffix = {'raw_stim','raw_ctrl','bs_stim','bs_ctrl'};
%         mkdir(fullfile(save_dir, 'avg_traces'));
%         saveas(fig_idx, fullfile(save_dir, 'avg_traces', strcat('avg_traces_', fig_suffix{fig_idx}, '.svg')));
%         saveas(fig_idx, fullfile(save_dir, 'avg_traces', strcat('avg_traces_', fig_suffix{fig_idx}, '.fig')));
%     end
% end
