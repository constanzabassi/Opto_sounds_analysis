function all_trial_info = get_trial_info_contexts(info,context_type,savepath)
%GOAL: get virmen trial info for each trial relative to bad_frames!
%old function took imaging_st{1,dataset_index} and alignment frames from
%opto_alignfunction

%now will just reload it for every dataset
%
    all_trial_info = struct('mouse_date', {}, 'serverid', {}, 'opto', {}, 'ctrl', {});

    for dataset_index = 1:length(info.mouse_date)
        disp(['Processing mouse ', num2str(dataset_index)]);
        
        % Retrieve mouse and server ID
        current_dataset = info.mouse_date{dataset_index};
        current_server = info.serverid{dataset_index};

        % Load bad frames
%         load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'bad_frames.mat'));
        load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'context_stim', '60', 'context_tr.mat')); %trials separated by context
        load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'context_stim', '60', 'bad_frames.mat'));
%         load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'context_stim', '60', 'exp.mat'));
%         load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'context_stim', '60', 'nonexp.mat'));
        load('V:\Connie\results\opto_sound_2025\context\data_info\alignment_frames.mat'); %using bad frames exp and non exp to get these from opto_align function (from pooled_data sound_opto code)

        %load alignment data
        load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'alignment_info.mat'));
        %%% 2) Select imaging and bad_frames!
        % Depending on context_type, the loaded variable name may differ.
        if strcmpi(context_type, 'passive')
            load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'passive','imaging.mat'));
        else
            load(fullfile(num2str(current_server), 'Connie', 'ProcessedData', num2str(current_dataset), 'VR','imaging.mat'));
        end

        %get frame lengths across t series files to make sure frame indices
        %make sense!
        frame_lengths = [];
        frame_lengths = cellfun(@length,{alignment_info.frame_times}); %across all imaged files 
        frame_lengths = [0,cumsum(frame_lengths)];

        %   get valid trials with imaging data
%         imaging = imaging_st{1, dataset_index};
        empty_trials = find(cellfun(@isempty, {imaging.good_trial}));
        good_trials = setdiff(1:length(imaging), empty_trials); % Only trials with all imaging data considered
        if strcmpi(context_type, 'passive')
            [~,alignment_frames_single_mouse,~,~] = find_align_info (imaging,30,2); %gets alignment frames for trial events (stimulus, turn, iti etc)
        else
            [~,alignment_frames_single_mouse,~,~] = find_align_info (imaging,30);
        end

        % Initialize trial info containers
        trial_info_opto = {};
        trial_info_ctrl = {};
        count = 0;
        count2 = 0;

        % Process opto trials
        for trial = good_trials
            count2 = count2+1;
            frames_to_add = alignment_frames_single_mouse(1,count2)+frame_lengths(imaging(trial).file_num)-1; %first row is sound onset which is opto onset!
            if strcmpi(context_type, 'passive')
                frames_to_add = alignment_frames_single_mouse(1,count2);
            end
            matched_trial = find(abs(imaging(trial).frame_id(1)+frames_to_add - alignment_frames{1, dataset_index}.opto(:, 61)) < 4); % Onset of opto has to be within 4 frames max!
            if ~isempty(matched_trial)
                count = count+1;
                %trial_info_opto{end + 1} = imaging(trial).virmen_trial_info;
                trial_info_opto(count).correct = imaging(trial).virmen_trial_info.correct;
                trial_info_opto(count).left_turn = imaging(trial).virmen_trial_info.left_turn;
                trial_info_opto(count).condition = imaging(trial).virmen_trial_info.condition ;
                trial_info_opto(count).is_stim_trial = imaging(trial).virmen_trial_info.is_stim_trial ;
                trial_info_opto(count).trial_id = trial; %relative to virmen
                trial_info_opto(count).matched_id = matched_trial; %relative to opto trials
                trial_info_opto(count).frame_onset = imaging(trial).frame_id(1)+frames_to_add;
            end
        end

        % Process control trials
        count = 0;
        count2 = 0;
        % Process opto trials
        for trial = good_trials
            count2 = count2+1;
            frames_to_add = alignment_frames_single_mouse(1,count2)+frame_lengths(imaging(trial).file_num)-1; %first row is sound onset which is opto onset!
            if strcmpi(context_type, 'passive')
                frames_to_add = alignment_frames_single_mouse(1,count2);
            end
            matched_trial = find(abs(imaging(trial).frame_id(1)+frames_to_add - alignment_frames{1, dataset_index}.ctrl(:, 61)) < 4); % Onset of opto
            if ~isempty(matched_trial)
                count = count+1;
                %trial_info_opto{end + 1} = imaging(trial).virmen_trial_info;
                trial_info_ctrl(count).correct = imaging(trial).virmen_trial_info.correct;
                trial_info_ctrl(count).left_turn = imaging(trial).virmen_trial_info.left_turn;
                trial_info_ctrl(count).condition = imaging(trial).virmen_trial_info.condition ;
                trial_info_ctrl(count).is_stim_trial = imaging(trial).virmen_trial_info.is_stim_trial ;
                trial_info_ctrl(count).trial_id = trial;
                trial_info_ctrl(count).matched_id = matched_trial;
                trial_info_ctrl(count).frame_onset = imaging(trial).frame_id(1)+frames_to_add;
            end
        end
        
        % Store the collected trial info into all_trial_info
        all_trial_info(dataset_index).mouse_date = current_dataset;
        all_trial_info(dataset_index).serverid = current_server;
        all_trial_info(dataset_index).opto = trial_info_opto;
        all_trial_info(dataset_index).ctrl = trial_info_ctrl;
    end

    %save the output?
if ~isempty(savepath)
    outdir = fullfile(savepath);
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    save(fullfile(outdir,[context_type '_all_trial_info.mat']), 'all_trial_info');
end
end