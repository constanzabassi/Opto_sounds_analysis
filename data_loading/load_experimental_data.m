function [data,context_tr] = load_experimental_data(server, dataset_path, path_string, multiple_contexts)
%load neural data (dff/deconv), stimulus onsets, and context_trials if
%multiple trials available

    base_path = fullfile(server, 'Connie', 'ProcessedData', dataset_path);
    
    % Load core files
    data.bad_frames = load_mat_file(base_path, 'bad_frames.mat');
    data.deconv = load_mat_file(fullfile(base_path, 'deconv'), 'deconv.mat');
    data.dff = load_mat_file(base_path, 'dff.mat');
    
    % Load experiment files
    exp_path = fullfile(base_path, path_string);
    data.exp = load_mat_file(exp_path, 'exp.mat'); %opto trials
    data.nonexp = load_mat_file(exp_path, 'nonexp.mat'); %control trials!

    % Load context data if required
    if multiple_contexts
        %has the trials relative to bad frames and stim/ctrl for each context
        context_tr = load_mat_file(exp_path, 'context_tr.mat');
    end
end

