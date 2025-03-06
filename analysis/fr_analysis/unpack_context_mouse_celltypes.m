function [deconv_response,chosen_cells] = unpack_context_mouse_celltypes(deconv_st,sig_mod_boot,all_celltypes,varagin)
% figure(96);clf
% t = tiledlayout(1,size(deconv_st,1),'TileSpacing','Compact','Padding','Compact');
fieldss = fields(all_celltypes{1,1});

if nargin > 3
    chosen_mice = varagin;
else
    chosen_mice = 1:size(deconv_st,2);
end

for mouse = chosen_mice
    chosen_cells ={};
    if ~isempty(sig_mod_boot) && ~isempty(sig_mod_boot{1,mouse})
        for ce = 1:length(fieldss)
            chosen_cells{ce} = sig_mod_boot{1,mouse}(find(ismember(sig_mod_boot{1,mouse},all_celltypes{1,mouse}.(fieldss{ce})))) ;
%             chosen_cells{2} =  sig_mod_boot{1,mouse}(find(ismember(sig_mod_boot{1,mouse},all_celltypes{1,mouse}.pv_cells)));
%             chosen_cells{3} =  sig_mod_boot{1,mouse}(find(ismember(sig_mod_boot{1,mouse},all_celltypes{1,mouse}.som_cells)));
        end
    elseif isempty(sig_mod_boot)
        for ce = 1:length(fieldss)
            chosen_cells{ce} = all_celltypes{1,mouse}.(fieldss{ce});
%             chosen_cells{2} = all_celltypes{1,mouse}.pv_cells;
%             chosen_cells{3} = all_celltypes{1,mouse}.som_cells;
        end
    end
    
    for context = 1:size(deconv_st,1)
        for cel = 1:length(fieldss)
            if ~isempty(chosen_cells) && length(chosen_cells{cel})>1 && all(cellfun(@(x) size(x.stim,1), {deconv_st{1:size(deconv_st,1),mouse}}) > 2) %at least 2 cells of this cell type!, at least 3 trials across all contexts for this mouse! % 
                deconv_response{context,mouse,cel}.stim = deconv_st{context,mouse}.stim(:,chosen_cells{cel},:);
                deconv_response{context,mouse,cel}.ctrl = deconv_st{context,mouse}.ctrl(:,chosen_cells{cel},:);
            else
                deconv_response{context,mouse,cel}.stim = nan;
                deconv_response{context,mouse,cel}.ctrl = nan;
            end
        end
    end
end


%     nexttile
%     title(celltypes_ids{2,cel})
%     for cel = 1:length(celltypes_ids)
%         for context = 1:size(deconv_st,1)
%             hold on
%             a(context) = plot(binss,stim_cdf);
%             set( a(context), 'LineWidth', 1.5, 'LineStyle', lineStyles_contexts{context}, 'color',colors(cel,:));%linecolors(2,:));
%         
%             b(context) = plot(binss,ctrl_cdf);
%             set( b(context), 'LineWidth', 1.5, 'LineStyle', lineStyles_contexts{context}, 'color',[0.5 0.5 0.5]);%linecolors(2,:));
%     
%             grid on
%     
%         
%         end
%     end
