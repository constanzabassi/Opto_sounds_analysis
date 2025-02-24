function params  = experiment_config()
        % Configuration settings for optical imaging analysis
    % Get base experiment parameters
    params.info = get_info_params();
    params.frames = get_frame_params();
    params.mod = get_modulation_params();
    params.mod_sounds = get_modulation_params_sounds()
    params.selectivity_sounds = get_selectivity_params_sounds()
end

function info = get_info_params()
    % Dataset and path information
    info.mouse_date = {'HA11-1R/2023-05-05'	'HA11-1R/2023-04-13'	'HA2-1L/2023-04-12'	'HA2-1L/2023-05-05'	'HA1-00/2023-06-29'	'HA1-00/2023-08-28'	'HE4-1L1R/2023-08-21'	'HE4-1L1R/2023-08-24'	'HA10-1L\2023-04-10'	'HA10-1L\2023-04-17'	'HA10-1L\2023-04-12'	'HA11-1R\2023-04-07'	'HA11-1R\2023-05-01'	'HA11-1R\2023-05-02'	'HA2-1L\2023-04-28'	'HA2-1L\2023-05-01'	'HA1-00\2023-06-27'	'HA1-00\2023-07-07'	'HA1-00\2023-08-25'	'HE4-1L1R\2023-08-14'	'HE4-1L1R\2023-08-28'	'HE4-1L1R\2023-09-04'	'HE4-1L1R\2023-09-11'	'HA10-1L\2023-03-31' 'HE1-00\2023-05-30'};
    info.serverid = {'V:'	'V:'	'V:'	'V:'	'V:'	'W:'	'W:'	'W:'	'V:'	'V:'	'V:'	'V:'	'V:'	'V:'	'V:'	'V:'	'V:'	'W:'	'W:'	'W:'	'W:'	'W:'	'W:'	'V:'    'V:'};
    info.path_string = 'context_stim/updated'; %used get_vaid_stim_onsets_trials to get trials with imaging structure info (from context_stim/60)
    info.savepath = 'V:/Connie/results/opto_sound_2025/context';
    info.savepath_sounds = 'V:/Connie/results/opto_sound_2025/context/sounds';
    mkdir(info.savepath);
    
    % Ensure 1xN format for cell arrays
    info.mouse_date = reshape(info.mouse_date, 1, []);
    info.serverid = reshape(info.serverid, 1, []);
end

function frames = get_frame_params()
    % Frame indices for analysis
    frames.before = 50:60;
    frames.after = 63:92;
    frames.stim = frames.before(end) + 1;
    frames.window = [60, 60];
end

function mod = get_modulation_params()
    % Modulation analysis parameters
    mod.mod_type = 'ctrl'; % 'prepost', 'ctrl', 'influence', 'prepost_sound'
    mod.mode = 'separate'; %'pooled', 'separate', 'simple'
    mod.response_range = {63:92; 51:60};
    mod.data_type_dff = 1; %dff or deconv
    mod.nShuffles = 10000; %total shuffles for bootstrapping
    mod.simple_or_not = 0; %if simple takes given trials and uses those to compute (no balancing or separating into left and right)
end

function mod = get_modulation_params_sounds()
    % Modulation analysis parameters
    mod.mod_type = 'prepost_sound'; % 'prepost', 'ctrl', 'influence', 'prepost_sound'
    mod.mode = 'separate'; %'pooled', 'separate', 'simple'
    mod.response_range = {63:92; 51:60};
    mod.data_type_dff = 1; %dff or deconv
    mod.nShuffles = 10000; %total shuffles for bootstrapping
    mod.simple_or_not = 0; %if simple takes given trials and uses those to compute (no balancing or separating into left and right)
end

function mod = get_selectivity_params_sounds()
    % Modulation analysis parameters
    mod.mod_type = 'ctrl'; % 'prepost', 'ctrl', 'influence', 'prepost_sound'
    mod.mode = 'selectivity'; %'pooled', 'separate', 'simple'
    mod.response_range = {63:92; 51:60};
    mod.data_type_dff = 1; %dff or deconv
    mod.nShuffles = 10000; %total shuffles for bootstrapping
    mod.simple_or_not = 0; %if simple takes given trials and uses those to compute (no balancing or separating into left and right)
end