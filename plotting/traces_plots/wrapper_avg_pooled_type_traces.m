function wrapper_avg_pooled_type_traces(context_data,all_celltypes,chosen_mice,savepath,data_type,plot_info,baseline)

%get distinct pools
sig_cells = [];

%get context,mouse,celltype responses (across all trials (not
%separated by left or rigth)- so overall avg)
[neural_response,~] = unpack_context_mouse_celltypes(context_data,sig_cells,all_celltypes,chosen_mice); %context_data.deconv_interp

%plot avg traces (plotting active and passive)
avg_across_neurons = 0; %SEM across all neurons vs across datasets
contexts_to_plot = [1,2]; %

%plot can include baseline subtraction but right now took it out
plot_info.trace_modes = {'raw', 'bs'}; %{'raw', 'bs'}
plot_avg_traces_baseline_subtracted(neural_response(contexts_to_plot,:,:),plot_info.colors_pooled_3contexts,{'-','-'},plot_info.pooled_names,1:122,[60,63],savepath,avg_across_neurons,[data_type num2str(baseline(1)) 'to' num2str(baseline(end))],plot_info,baseline);
