function [modl_fit, index_updated,index2_updated,stats] = cdf_index_sigcell(sig_mod_boot, all_celltypes, index, plot_info, save_path, string1, string2,index_lab, varargin)
index_updated = {};
index2_updated = {};
positions = utils.calculateFigurePositions(1, 5, .5, []);
    if nargin > 8
        minmax = varargin{1,1};
    else
        minmax = [-1, 1];
    end

    % Create figure
figure(2324); clf;

positions = utils.calculateFigurePositions(1,5,.5,[]);
mainAx = axes('Position',positions(1,:));
hold(mainAx,'on')

modl_fit = cell(1,3);
celltype_fields = fields(all_celltypes{1});
n_celltypes = length(celltype_fields);

for cell_type = 1:n_celltypes

    all_x = [];
    all_y = [];

    for dataset_index = 1:length(all_celltypes)

        if ~isempty(sig_mod_boot)
            selected_cells = sig_mod_boot{dataset_index}(ismember(sig_mod_boot{dataset_index}, ...
                all_celltypes{dataset_index}.(celltype_fields{cell_type})));
        else
            if size(all_celltypes,1) > 1
                selected_cells = all_celltypes{ctx,dataset_index}.(celltype_fields{cell_type});
            else
                selected_cells = all_celltypes{dataset_index}.(celltype_fields{cell_type});
            end
        end

        x = index{dataset_index,1}(selected_cells);
        y = index{dataset_index,2}(selected_cells);

        all_x = [all_x; x(:)];
        all_y = [all_y; y(:)];

    end

    if ~isempty(all_x)

        % compute CDF bins
        binss = linspace(minmax(1),minmax(2),25);

        % CDF X
        [cdf_x,~] = make_cdf(all_x,binss);

        % CDF Y
        [cdf_y,~] = make_cdf(all_y,binss);

        % plot
        plot(mainAx,binss,cdf_x,...
            'Color',plot_info.colors_celltypes(cell_type,:),...
            'LineWidth',1)

        plot(mainAx,binss,cdf_y,...
            'Color',plot_info.colors_celltypes(cell_type,:),...
            'LineStyle',':',...
            'LineWidth',1)

        % stats
        [r,p_val] = corr(all_x,all_y);
         field_name = regexprep(string1, ' ', '');
            field_name2 = regexprep(string2, ' ', '');
            field_name = regexprep(field_name, '[^a-zA-Z0-9]', '_');
            field_name2 = regexprep(field_name2, '[^a-zA-Z0-9]', '_');
            field_name = regexprep(field_name, '_+', '_');
            field_name2 = regexprep(field_name2, '_+', '_');
            field_name = regexprep(field_name,  '^_|_$', '');
            field_name2 = regexprep(field_name2,  '^_|_$', '');
            field_name = regexprep(field_name,  'Δ', '');
            field_name2 = regexprep(field_name2,  'Δ', '');
            
            
            % Ensure it doesn’t start with a number
            if ~isletter(field_name(1))
                field_name = ['f_' field_name];
            end
            if ~isletter(field_name2(1))
                field_name2 = ['f_' field_name2];
            end
            all_field_names = strcat(field_name,field_name2);

        stats.(celltype_fields{cell_type}).(field_name) = get_basic_stats(all_x);
        stats.(celltype_fields{cell_type}).(field_name2) = get_basic_stats(all_y);

        stats.(celltype_fields{cell_type}).corr.r = r;
        stats.(celltype_fields{cell_type}).corr.p_val = p_val;

        index_updated{cell_type} = all_x';
        index2_updated{cell_type} = all_y';

    end
end

xlabel(mainAx,index_lab)
ylabel(mainAx,'CDF')

xlim([minmax(1) minmax(2)])
ylim([0 1])

box off
set(mainAx,'FontSize',7)

% legend({string1,string2})
    set(mainAx, 'FontSize', 7)

    % Use painters renderer (vector-safe)
    set(gcf, 'Renderer', 'painters');
    
%     set(gcf, 'Position', [100, 100, 200, 200]);  % [left bottom width height]
    set(mainAx, 'FontSize', 7, 'Units', 'inches', 'Position', positions(1, :));
    mainAx.XLabel.FontSize = 7;
    mainAx.YLabel.FontSize = 7;

   xline(0,  '--', 'Color', [0.5 0.5 0.5]);
        


    % Save figure
    sig_cel_string = ~isempty(sig_mod_boot);
    if ~isempty(save_path)
        if ~exist(save_path, 'dir')
            mkdir(save_path);
        end
        string1 = strrep(string1, '\Delta', 'Δ');
        string2 = strrep(string2, '\Delta', 'Δ');

        string2 = strrep(string2, '\', '');
        string1 = strrep(string1, '\', '');
        saveas(gcf, fullfile(save_path, ['cdf_index_sigcells_histogram' num2str(sig_cel_string) '_' string1 '_' string2 '.png']));
        saveas(gcf, fullfile(save_path, ['cdf_index_sigcells_histogram' num2str(sig_cel_string) '_' string1 '_' string2 '.svg']));
        exportgraphics(gcf, fullfile(save_path, ['cdf_index_sigcells_histogram' num2str(sig_cel_string) '_' string1 '_' string2 '.pdf']), 'ContentType', 'vector');
        %save_stats!
        save(fullfile(save_path,['stats_cdf_index_sigcells_histogram' num2str(sig_cel_string) '_' string1 '_' string2 '.mat']),'stats');





    end
end
