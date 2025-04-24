
function plot_avg_heatmap_direction_comparison(avg_results, selectivity_results,save_dir,varargin)
    % Plot average responses for both left and right selective populations
    
    
        params = get_default_params();

    
    % Create figure with 2 rows (left/right selective) x 2 columns (left/right sounds)
    figure('Position', [100 100 900 800]);

        params.colors.active = [0,0,0];
    params.colors.passive = [0.5,0.5,0.5];
    params.stim_onset = 61;

    if nargin > 3
        % Plot left selective population
        save_name = plot_population('both.left', 1, varargin);
        
        % Plot right selective population
        save_name = plot_population('both.right', 2, varargin);

        % Plot non-selective population (added as the third row)
        save_name = plot_population('both.nonsel', 3, varargin);
    else
        % Plot left selective population
        save_name = plot_population('both.left', 1);
        
        % Plot right selective population
        save_name = plot_population('both.right', 2);

        % Plot non-selective population (added as the third row)
        save_name = plot_population('both.nonsel', 3);
    end

    
    sgtitle('Population Responses by Selectivity', 'FontWeight', 'normal');
    
    function save_name = plot_population(pool_type, row,varargin)
        parts = strsplit(pool_type, '.');
        if strcmp(parts{2}, 'nonsel')
            cell_indices = selectivity_results.both.nonsel.cell_indices;
        else
            cell_indices = selectivity_results.(parts{1}).(parts{2}).cell_indices;
        end
        directions = {'Left', 'Right'};
        
        for dir = 1:2
            subplot(3, 2, (row-1)*2 + dir);
            hold on;
            
            % Get responses for this direction in both contexts
            active_responses = avg_results{1,1}.(lower(directions{dir})).neuron_mean(cell_indices,:);
            passive_responses = avg_results{1,2}.(lower(directions{dir})).neuron_mean(cell_indices,:);
                        
            % Plot with shadedErrorBar
            time_vector = 1:size(active_responses,2);
            mid = ceil(size([active_responses;passive_responses],1)*0.05);
            imagesc(([active_responses;zeros(mid,size(active_responses,2));passive_responses]));

            
            % Add stimulus onset line
            xline(params.stim_onset, '-w', 'LineWidth', 1);
            utils.set_current_fig
            
            % Format plot
            xlabel('Time (s)');
            if nargin > 2 && numel(varargin)>0 &&length(varargin{1, 1}{1, 1}{1, 1})>1
                ylabel(varargin{1, 1}{1, 1}{1, 1}{1, 2});
            else
                ylabel('Response Nuerons (ΔF/F)');
            end
                        % Sanitize ylabel text for filenames (remove spaces and special characters)

            title(sprintf('%s Sound (n=%d)', directions{dir}, length(cell_indices)));
            
            % Add axis ticks
            [xticks_in, xticks_lab] = utils.x_axis_sec_aligned(...
                params.stim_onset, length(time_vector), 1);
            xticks(xticks_in);
            xticklabels(xticks_lab);
            caxis([-.02,.4])
                        
%             if nargin > 2 && numel(varargin)>0
%                 ylim(varargin{1, 1}{1, 1}{1, 1}{1});
%             else
%                 ylim([.10 .4])
%             end
        end

        % Add row label
        if strcmp(parts{2}, 'nonsel')
            row_label = 'Non-Selective';
        else
            row_label = sprintf('%s Selective', parts{2});
        end
        
        % Add row label
        annotation('textbox', [0.01, 1.05 - row * 0.33, 0.1, 0.1], ...
            'String', row_label, ...
            'EdgeColor', 'none', 'FontSize', 12);

        save_name = 'test';
    end

    
    % Save plots with the ylabel in the filename
    if ~isempty(save_dir)
        mkdir(save_dir);
        saveas(gcf, fullfile(save_dir, ['avg_heatmap_direction_comparison_' save_name '.png']));
        saveas(gcf, fullfile(save_dir, ['avg_heatmap_direction_comparison_' save_name '.fig']));
        exportgraphics(gcf, fullfile(save_dir, ['avg_heatmap_direction_comparison_' save_name '.pdf']), 'ContentType', 'vector');
    end
end

function params = get_default_params()
    params = struct();
    params.colors.active = [0.86, 0.15, 0.49];   % Pink
    params.colors.passive = [0.39, 0.56, 1.0];   % Light blue
    params.patch_saturation = 0.2;
    params.stim_onset = 61;
end