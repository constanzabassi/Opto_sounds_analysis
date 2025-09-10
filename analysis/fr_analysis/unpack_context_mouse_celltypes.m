function [deconv_response,chosen_cells] = unpack_context_mouse_celltypes(deconv_st,sig_mod_boot,all_celltypes,varagin)
%reconfigure structure into cell array of
%neural_data{context,mouse,cel}.stim or neural_data{context,mouse,cel}.ctrl

fieldss = fields(all_celltypes{1,1});

if nargin > 3
    chosen_mice = varagin;
else
    chosen_mice = 1:size(deconv_st,2);
end
chosen_cells ={};
for dataset_index = chosen_mice
    dataset_index
    if dataset_index == 3
        a=1;
    end
    if ~isempty(sig_mod_boot) && ~isempty(sig_mod_boot{1,dataset_index})
        for ce = 1:length(fieldss)
            chosen_cells{dataset_index,ce} = sig_mod_boot{1,dataset_index}(find(ismember(sig_mod_boot{1,dataset_index},all_celltypes{1,dataset_index}.(fieldss{ce})))) ;
        end
    elseif isempty(sig_mod_boot)
        for ce = 1:length(fieldss)
            chosen_cells{dataset_index,ce} = all_celltypes{1,dataset_index}.(fieldss{ce});
        end
    else isempty(sig_mod_boot{1,dataset_index})
        for ce = 1:length(fieldss)
        chosen_cells{dataset_index,ce} = {};
        end
    end
    
    for context = 1:size(deconv_st,1)
        context
        for cel = 1:length(fieldss)
            if ~isempty(chosen_cells) && length(chosen_cells{dataset_index,cel})>1 && all(cellfun(@(x) size(x.stim,1), {deconv_st{1:size(deconv_st,1),dataset_index}}) > 2) %at least 2 cells of this cell type!, at least 3 trials across all contexts for this mouse! % 
                deconv_response{context,dataset_index,cel}.stim = deconv_st{context,dataset_index}.stim(:,chosen_cells{dataset_index,cel},:);
                deconv_response{context,dataset_index,cel}.ctrl = deconv_st{context,dataset_index}.ctrl(:,chosen_cells{dataset_index,cel},:);
                if size(all_celltypes,1)>1 && ~isempty(all_celltypes{context,dataset_index}.(fieldss{cel}))
                    chosen_cells{dataset_index,cel} = all_celltypes{context,dataset_index}.(fieldss{cel});
                    deconv_response{context,dataset_index,cel}.stim = deconv_st{context,dataset_index}.stim(:,chosen_cells{dataset_index,cel},:);
                    deconv_response{context,dataset_index,cel}.ctrl = deconv_st{context,dataset_index}.ctrl(:,chosen_cells{dataset_index,cel},:);
                end
            else
                deconv_response{context,dataset_index,cel}.stim = nan;
                deconv_response{context,dataset_index,cel}.ctrl = nan;
            end
        end
    end
end