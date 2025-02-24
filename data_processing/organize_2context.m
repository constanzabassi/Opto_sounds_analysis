function [dff_context,stim_trials_context,ctrl_trials_context] = organize_2context(contex_st,context2_st,varargin)
for m = 1:length(contex_st)
    condition1_trials_all = [];
    condition2_trials_all = [];
    for c = 1:nargin
        if c == 1
            stim_trials = 1:size(contex_st{1,m}.stim,1);
            dff_context{c,m}.stim = contex_st{1,m}.stim;
            ctrl_trials = 1:size(contex_st{1,m}.ctrl,1);
            dff_context{c,m}.ctrl = contex_st{1,m}.ctrl;
            condition1_trials_all = [condition1_trials_all,{stim_trials}];
            condition2_trials_all = [condition2_trials_all,{ctrl_trials}];
        else
            stim_trials = stim_trials(end)+1:stim_trials(end)+size(context2_st{1,m}.stim,1);
            dff_context{c,m}.stim = context2_st{1,m}.stim;
            ctrl_trials = ctrl_trials(end)+1:ctrl_trials(end)+size(context2_st{1,m}.ctrl,1);
            dff_context{c,m}.ctrl = context2_st{1,m}.ctrl;
            condition1_trials_all = [condition1_trials_all,{stim_trials}];
            condition2_trials_all = [condition2_trials_all,{ctrl_trials}];
        end
    end
    stim_trials_context{m} = condition1_trials_all;
    ctrl_trials_context{m} = condition2_trials_all;
end