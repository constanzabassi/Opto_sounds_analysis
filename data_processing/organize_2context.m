function [dff_context, stim_trials_context, ctrl_trials_context] = organize_2context(context_st, context2_st, varargin)
% ORGANIZE_2CONTEXT organizes neural activity data for two different contexts.
% It extracts and structures stimulus (stim) and control (ctrl) trial data 
% from two context structures and ensures trial indices are managed properly.
%
% Inputs:
%   context_st  - Cell array containing neural data for the first context.
%   context2_st - Cell array containing neural data for the second context.
%   varargin    - Optional additional contexts (not currently used but allows extension).
%
% Outputs:
%   dff_context         - Structure array with organized neural data.
%   stim_trials_context - Cell array of stimulus trial indices for each dataset.
%   ctrl_trials_context - Cell array of control trial indices for each dataset.
%
% Author: CB 03/04/2025

for current_dataset = 1:length(context_st)
    condition1_trials_all = [];
    condition2_trials_all = [];
    
    for current_context = 1:nargin
        if current_context == 1  % First context
            % Get trial indices
            stim_trials = 1:size(context_st{1, current_dataset}.stim, 1);
            ctrl_trials = 1:size(context_st{1, current_dataset}.ctrl, 1);
            
            % Store data in output structure
            dff_context{current_context, current_dataset}.stim = context_st{1, current_dataset}.stim;
            dff_context{current_context, current_dataset}.ctrl = context_st{1, current_dataset}.ctrl;
            
            % Store trial indices
            condition1_trials_all = [condition1_trials_all, {stim_trials}];
            condition2_trials_all = [condition2_trials_all, {ctrl_trials}];
            
            % Check if z-scored data exists and store it
            if isfield(context_st{1, current_dataset}, 'z_stim') && isfield(context_st{1, current_dataset}, 'z_ctrl')
                dff_context{current_context, current_dataset}.z_stim = context_st{1, current_dataset}.z_stim;
                dff_context{current_context, current_dataset}.z_ctrl = context_st{1, current_dataset}.z_ctrl;
            end
        else  % Second context (or additional contexts)
            % Adjust trial indices to avoid overlap
            stim_trials = stim_trials(end) + (1:size(context2_st{1, current_dataset}.stim, 1));
            ctrl_trials = ctrl_trials(end) + (1:size(context2_st{1, current_dataset}.ctrl, 1));
            
            % Store data in output structure
            dff_context{current_context, current_dataset}.stim = context2_st{1, current_dataset}.stim;
            dff_context{current_context, current_dataset}.ctrl = context2_st{1, current_dataset}.ctrl;
            
            % Store trial indices
            condition1_trials_all = [condition1_trials_all, {stim_trials}];
            condition2_trials_all = [condition2_trials_all, {ctrl_trials}];
            
            % Check if z-scored data exists and store it
            if isfield(context2_st{1, current_dataset}, 'z_stim') && isfield(context2_st{1, current_dataset}, 'z_ctrl')
                dff_context{current_context, current_dataset}.z_stim = context2_st{1, current_dataset}.z_stim;
                dff_context{current_context, current_dataset}.z_ctrl = context2_st{1, current_dataset}.z_ctrl;
            end
        end
    end
    
    % Store trial indices for both conditions
    stim_trials_context{current_dataset} = condition1_trials_all;
    ctrl_trials_context{current_dataset} = condition2_trials_all;
end
