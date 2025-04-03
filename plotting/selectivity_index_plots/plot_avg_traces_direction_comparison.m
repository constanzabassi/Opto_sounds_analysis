
function plot_avg_traces_direction_comparison(avg_results, selectivity_results,save_dir,varargin)
    % Plot average responses for both left and right selective populations
    
    
        params = get_default_params();

    
    % Create figure with 2 rows (left/right selective) x 2 columns (left/right sounds)
    figure('Position', [100 100 900 800]);

        params.colors.active = 'k';
    params.colors.passive = [0.5,0.5,0.5];
    params.stim_onset = 61;

    if nargin > 3
        % Plot left selective population
        save_name = plot_population('both.left', 1, varargin);
        
        % Plot right selective population
        save_name = plot_population('both.right', 2, varargin);
    else
        % Plot left selective population
        save_name = plot_population('both.left', 1);
        
        % Plot right selective population
        save_name = plot_population('both.right', 2);
    end

    
    sgtitle('Population Responses by Selectivity', 'FontWeight', 'normal');
    
    function save_name = plot_population(pool_type, row,varargin)
        parts = strsplit(pool_type, '.');
        cell_indices = selectivity_results.(parts{1}).(parts{2}).cell_indices;
        directions = {'Left', 'Right'};
        
        for dir = 1:2
            subplot(2, 2, (row-1)*2 + dir);
            hold on;
            
            % Get responses for this direction in both contexts
            active_responses = avg_results{1,1}.(lower(directions{dir})).neuron_mean(cell_indices,:);
            passive_responses = avg_results{1,2}.(lower(directions{dir})).neuron_mean(cell_indices,:);
            
            % Calculate mean and SEM
            active_mean = mean(active_responses, 1);
            passive_mean = mean(passive_responses, 1);
            active_sem = std(active_responses, [], 1)/sqrt(size(active_responses,1));
            passive_sem = std(passive_responses, [], 1)/sqrt(size(passive_responses,1));
            
            % Plot with shadedErrorBar
            time_vector = 1:size(active_responses,2);
            h1 = shadedErrorBar(time_vector, active_mean, active_sem, ...
                'lineProps', {'Color', params.colors.active, 'LineWidth', 1.5}, ...
                'patchSaturation', 0.5);
            h2 = shadedErrorBar(time_vector, passive_mean, passive_sem, ...
                'lineProps', {'Color', params.colors.passive, 'LineWidth', 1.5}, ...
                'patchSaturation', 0.5);
            
            % Add stimulus onset line
            xline(params.stim_onset, '--k', 'LineWidth', 1);
            
            % Format plot
            xlabel('Time (s)');
            if nargin > 2 && numel(varargin)>0 &&length(varargin{1, 1}{1, 1}{1, 1})>1
                ylabel(varargin{1, 1}{1, 1}{1, 1}{1, 2});
                ylabel_text = varargin{1, 1}{1, 1}{1, 1}{1, 2};
            else
                ylabel('Response (ΔF/F)');
                ylabel_text = 'Response (ΔF/F)';
            end
                        % Sanitize ylabel text for filenames (remove spaces and special characters)
            save_name = strrep(ylabel_text, ' ', '_');  % Replace spaces with underscores
            save_name = regexprep(save_name, '[^\w_]', '');  % Remove non-alphanumeric characters except underscores


            title(sprintf('%s Sound (n=%d)', directions{dir}, length(cell_indices)));
            
            % Add axis ticks
            [xticks_in, xticks_lab] = utils.x_axis_sec_aligned(...
                params.stim_onset, length(time_vector), 1);
            xticks(xticks_in);
            xticklabels(xticks_lab);
            
            % Add legend only once
            if row == 1 && dir == 1
                legend([h1.mainLine, h2.mainLine], {'Active', 'Passive'}, ...
                    'Location', 'best', 'box', 'off');
            end
            
            utils.set_current_fig;
            
            if nargin > 2 && numel(varargin)>0
                ylim(varargin{1, 1}{1, 1}{1, 1}{1});
            else
                ylim([.10 .4])
            end
        end
        
        % Add row label
        annotation('textbox', [0.01, 1.05-row*0.5, 0.1, 0.1], ...
            'String', sprintf('%s Selective', parts{2}), ...
            'EdgeColor', 'none', 'FontSize', 12);

        
    end

    % Save plots with the ylabel in the filename
    if ~isempty(save_dir)
        mkdir(save_dir);
        saveas(gcf, fullfile(save_dir, ['avg_traces_direction_comparison_' save_name '.png']));
        saveas(gcf, fullfile(save_dir, ['avg_traces_direction_comparison_' save_name '.fig']));
        exportgraphics(gcf, fullfile(save_dir, ['avg_traces_direction_comparison_' save_name '.pdf']), 'ContentType', 'vector');
    end
end

function params = get_default_params()
    params = struct();
    params.colors.active = [0.86, 0.15, 0.49];   % Pink
    params.colors.passive = [0.39, 0.56, 1.0];   % Light blue
    params.patch_saturation = 0.2;
    params.stim_onset = 61;
end