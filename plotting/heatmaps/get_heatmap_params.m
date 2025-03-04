function params = get_heatmap_params()
    params = struct();
    params.sort_method = 'peak';
    params.response_window = 62:122%92;  % 1s after stim onset
    params.clim = [-.5 1];           % Typical range for z-scores
    params.stim_onset = 61;
    params.colormap = 'viridis';
    params.sort_method = 'peak'; %'latency'/'peak'
    params.context_labels = {'Active','Passive','Spont'};
end