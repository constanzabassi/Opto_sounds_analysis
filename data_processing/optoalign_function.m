%addpath '\\136.142.49.216\runyan2\Connie\Code'
function [allcells,allcells_nogap,frames] = optoalign_function(stim, ctrl, bad_frames,dff,deconv, before_frames,after_frames)
z_dff=zscore(dff,0,2);
%number of frames before and after stim
x1=before_frames+1;
x2=1;
y1=2;
y2=after_frames+2;

intervals = bad_frames;
% intervals(:,1)= bad_frames(:,1)-1;
% intervals(:,2)= bad_frames(:,2)+2;

stim_matrix = [];
ctrl_matrix = [];
cellCount=size(dff,1);

stim_trials = [ctrl, stim]; %always ordered so control are first then exp!!! (not ordered but numerical values)
%keep stim_trials intervals data deconv nonexp exp before_frames after_frames
lenInt= length(stim_trials);%length includes only intervals in nonexp and exp
trials_stim=0;
trials_ctrl=0;
for cel=1:cellCount
    for i=1:lenInt 
       j = stim_trials(i);
       x= intervals(j,1);
       y= intervals(j,2);
       bfint= x-x1:x-x2; 
       afint= y+y1:y+y2;
       bfint2 = x+1:x+x1+1;%y-y2+1:y-1;
           if ismember( j, stim)
               trials_stim=trials_stim+1;
                allcells(cel).opto(trials_stim,:)= [dff(cel, [bfint afint])];
                allcells(cel).opto_deconv(trials_stim,:)= [deconv(cel, [bfint afint])];
                allcells(cel).z_opto(trials_stim,:)= [z_dff(cel, [bfint afint])];

                allcells_nogap(cel).opto(trials_stim,:)= [dff(cel,x-before_frames:x+(before_frames+1) )];%[bfint bfint2]%bfint2 afint
                allcells_nogap(cel).opto_deconv(trials_stim,:)= [deconv(cel, x-before_frames:x+(before_frames+1) )];
                allcells_nogap(cel).z_opto(trials_stim,:)= [z_dff(cel, x-before_frames:x+(before_frames+1) )];
                frames.opto(trials_stim,:) = [bfint afint];
           else 
                trials_ctrl=trials_ctrl+1;
                allcells(cel).control(trials_ctrl,:)= [dff(cel, [bfint afint])];
                allcells(cel).control_deconv(trials_ctrl,:)= [deconv(cel, [bfint afint])];
                allcells(cel).z_control(trials_ctrl,:)= [z_dff(cel, [bfint afint])];

                allcells_nogap(cel).control(trials_ctrl,:)= [dff(cel, x-before_frames:x+(before_frames+1) )];
                allcells_nogap(cel).control_deconv(trials_ctrl,:)= [deconv(cel, x-before_frames:x+(before_frames+1) )];
                allcells_nogap(cel).z_control(trials_ctrl,:)= [z_dff(cel, x-before_frames:x+(before_frames+1) )];
                frames.ctrl(trials_ctrl,:) = [bfint afint];
           end
    end
    trials_stim=0;
    trials_ctrl=0;
end

end %%