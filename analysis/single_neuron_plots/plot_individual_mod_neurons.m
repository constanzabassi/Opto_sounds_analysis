function plot_individual_mod_neurons(stim_data_subset, ctrl_data_subset, mod_index, sig_neurons, time_vector, dataset_index,left_trials_total,save_dir)
% PLOT_INDIVIDUAL_MOD_NEURONS generates plots of neural activity for selected neurons.
%
%   Inputs:
%       stim_data_subset  - 3D array of stimulation data [nTrials x nNeurons x nFrames].
%       ctrl_data_subset  - 3D array of control data [nTrials x nNeurons x nFrames] (optional).
%                           (In this example, we use only stimulation data.)
%       mod_index         - 1 x nNeurons vector of modulation indices.
%       sig_neurons       - Vector of neuron indices (local to the dataset) to plot.
%       time_vector       - Vector of time points corresponding to frames.
%                           If empty, frame indices are used.
%       save_dir          - (Optional) Directory where figures should be saved.
%
%   For each neuron in sig_neurons, the function creates a figure with:
%       - Subplot (2,1,1): an imagesc plot of the neuron's trial-by-trial activity.
%       - Subplot (2,1,2): a line plot of the mean response (averaged across trials).
%   The figure title includes the neuron ID and its modulation index.
%
%   Example:
%       plot_individual_mod_neurons(stim_data, ctrl_data, mod_index, sig_neurons, 1:120, 'C:\save_figs');
%
%   Author: Your Name, Date
    % If time_vector is empty, default to frame indices.
    if nargin < 5 || isempty(time_vector)
        time_vector = 1:size(stim_data_subset, 3);
    end
    if nargin < 6
        save_dir = '';
    end

    if length(sig_neurons) > 5
        total_to_plot = 5;
    else
        total_to_plot = length(sig_neurons);
    end
    % Loop over each selected neuron.
    for i = 1:total_to_plot
        neuron_id = sig_neurons(i);
        % Extract the neuron's data (across trials and frames).
        % This returns a matrix of size [nTrials x nFrames].
        neuron_data = squeeze(stim_data_subset(:, neuron_id, :));
        neuron_data_ctrl = squeeze(ctrl_data_subset(:, neuron_id, :));

        % Compute the mean response across trials.
        mean_response = mean(neuron_data, 1);
        mean_response_ctrl = mean(squeeze(ctrl_data_subset(:, neuron_id, :)));

        % Create a figure for this neuron.
        figure(i+(dataset_index(1)*10));clf;
        % Subplot 1: imagesc of trial-by-trial activity.
        subplot(2,1,1);
        imagesc([neuron_data;neuron_data_ctrl]); %stack trials
        colormap('viridis');  % Change colormap as desired.
%         colorbar;
        xlabel('Time');
        ylabel('Trial');
        xline(61,'-w','LineWidth',2)
        yline(left_trials_total(1),'-w','LineWidth',1) %size(neuron_data,1)/2
        yline(size(neuron_data,1),'--w','LineWidth',2)
        yline(size(neuron_data,1)+left_trials_total(2),'-w','LineWidth',1)
        [xticks_in,xticks_lab]  = utils.x_axis_sec_aligned(61,122,1);
        xticks(xticks_in);
        xticklabels(xticks_lab);

        caxis([-.3 1]);
        set_current_fig;
        set(gca,'FontSize',12);
        title(sprintf('Neuron %d Mod Index: %.3f', neuron_id, mod_index(neuron_id)),'FontWeight','normal');

        % Subplot 2: line plot of the mean response.
        subplot(2,1,2);
        hold on
        plot(time_vector, mean_response,'k', 'LineWidth', 1.5);
        plot(time_vector, mean_response_ctrl, 'LineWidth', 1.5,'Color', [0.5 0.5 0.5]);
        xlabel('Time');
        ylabel('Mean Response');
%         title('Mean Across Trials','FontWeight','normal');
        xlim([min(time_vector), max(time_vector)]);
        xline(61,'--k','LineWidth',1)
        [xticks_in,xticks_lab]  = utils.x_axis_sec_aligned(61,122,1);
        xticks(xticks_in);
        xticklabels(xticks_lab);
        set_current_fig;
        set(gca,'FontSize',12);
        hold off
        % Adjust figure size (optional).
        set(gcf, 'Position', [100, 100, 300, 400]);
        % Save the figure if a save directory is provided.
        if ~isempty(save_dir)
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            saveas(gcf, fullfile(save_dir, sprintf('Neuron_%d_dataset_%i_context_%c.fig', neuron_id,dataset_index(1),num2str(dataset_index(2)))));
            saveas(gcf, fullfile(save_dir, sprintf('Neuron_%d_dataset_%i_context_%c.png', neuron_id,dataset_index(1),num2str(dataset_index(2)))));
        end
    end
end