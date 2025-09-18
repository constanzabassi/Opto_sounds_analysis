function plot_linear_regression_lines(lme,tbl,context_all,ylabel_string,save_dir,varargin)
%extract table variables
var_names = tbl.Properties.VariableNames;

% separate into contexts for plotting
tbl_active = tbl(context_all == 0,var_names);
tbl_passive = tbl(context_all == 1,var_names);

% % Predict using model (used when context was a separate predictor
% [yhat_active,  yci_active]  = predict(lme, tbl_active);
% [yhat_passive, yci_passive] = predict(lme, tbl_passive);


[yhat_active,  yci_active]  = predict(lme, tbl);
tbl_active = tbl(:,var_names);


%% make plots!

positions = utils.calculateFigurePositions(1, 5, .5, []);

xvals = linspace(min(tbl{:,var_names{2}}), max(tbl{:,var_names{2}}), 100); %max value of the predictor axis (x)
% Passive line
pred_active = lme.Coefficients.Estimate(1) + ...
               lme.Coefficients.Estimate(2) * xvals;

figure(103); clf; hold on;
scatter(tbl{context_all==0,var_names{2}}, tbl{context_all==0,var_names{1}}, 5,'MarkerEdgeColor',[0.2 0.2 0.2], 'MarkerEdgeAlpha',.8)%[0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha',1)
scatter(tbl{context_all==1,var_names{2}}, tbl{context_all==1,var_names{1}}, 5,'MarkerEdgeColor',[0.8 0.8 0.8], 'MarkerEdgeAlpha',.8)%[0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha',1)

plot(xvals, pred_active, 'k', 'LineWidth', 2.2)


ylabel({ylabel_string;'(z-scored)'});

if nargin > 5
    xlabel(strcat({varargin{1,1},' Projection';'(z-scored)'}))
    ylabel_string_updated = [ylabel_string varargin{1,1}];
else
    xlabel({'Engagement Projection';'(z-scored)'})
end

%set figure
set(gca, 'FontSize', 7, 'Units', 'inches', 'Position', positions(1, :));
%include n, p value, r value, slope?
n_total = size(tbl,1);
[r,p_val] = corr(table2array(tbl(:,1)),table2array(tbl(:,2)));
utils.place_text_labels({['n = ', num2str(n_total)]},'k',0,5,'topleft',0.05)
utils.place_text_labels({['r = ', num2str(round(r,2))]},'k',0.1,5,'topleft',0.05)
utils.place_text_labels({['p = ',num2str(p_val, '%.1e')]},'k',0.2,5,'topleft',0.05)

utils.set_current_fig;

% SECOND PLOT WITHOUT ANY SCATTER POINTS
% Plot in figure 103
figure(1003); clf; hold on;
% Use actual predictor values (x) for each context
x_passive = tbl_passive{:, var_names{2}};
x_active = tbl_active{:, var_names{2}};

% Sort x and corresponding predictions for clean fill plots
[x_active_sorted, idx_active] = sort(x_active);
yhat_active_sorted = yhat_active(idx_active);
yci_active_sorted = yci_active(idx_active, :);

% Plot passive trial predictions and CIs
fill([x_active_sorted; flipud(x_active_sorted)], ...
     [yci_active_sorted(:,1); flipud(yci_active_sorted(:,2))], ...
     [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
plot(x_active_sorted, yhat_active_sorted, 'k','LineWidth', 2);

% Plot active trial predictions and CIs (if model includes context)

ylabel({ylabel_string;'(z-scored)'});
if nargin > 5
    xlabel(strcat({varargin{1,1},' Projection';'(z-scored)'}))
    ylabel_string_updated = [ylabel_string varargin{1,1}];
else
    ylabel_string_updated = ylabel_string;
    xlabel({'Engagement Projection';'(z-scored)'})
end
%set figure
set(gca, 'FontSize', 7, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

% Save results
if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(103,strcat('linear_regression_',num2str(ylabel_string_updated),'.fig'));
    exportgraphics(figure(103),strcat('linear_regression_',num2str(ylabel_string_updated),'.pdf'), 'ContentType', 'vector');

    saveas(1003,strcat('linear_lineonly_regression_',num2str(ylabel_string_updated),'.fig'));
    exportgraphics(figure(1003),strcat('linear_lineonly_regression_',num2str(ylabel_string_updated),'.pdf'), 'ContentType', 'vector');

end
