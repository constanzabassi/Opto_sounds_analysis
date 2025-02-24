function [stim_trials_combined, dff_combined] = combine_stim_info_dff_st(stim_info1, stim_info2,dff1, dff2)
stim_trials_combined = {};
for m = 1:length(stim_info1)
    stim_trials_combined{m,1} = [stim_info1{m,1};stim_info2{m,1}];
    stim_trials_combined{m,2} = [stim_info1{m,2},stim_info2{m,2}+length(stim_info1{m,1})];
    stim_trials_combined{m,3} = [stim_info1{m,3},stim_info2{m,3}+length(stim_info1{m,1})];

    dff_combined{1,m}.stim = [dff1{1,m}.stim; dff2{1,m}.stim];
    dff_combined{1,m}.ctrl = [dff1{1,m}.ctrl; dff2{1,m}.ctrl];
%     dff_combined{1,m}.z_stim = [dff1{1,m}.z_stim; dff2{1,m}.z_stim];
%     dff_combined{1,m}.z_ctrl = [dff1{1,m}.z_ctrl; dff2{1,m}.z_ctrl];

end