function plot_me_regression_lines(lme,engagement_proj_all,sound_proj_all,context_all,ylabel_string,save_dir,varargin)
% Create prediction table for both contexts
xvals = linspace(min(engagement_proj_all), max(engagement_proj_all), 100)';
response_var = strcat(ylabel_string, 'Proj');  % e.g., 'SoundProj'

% Ensure Context is categorical and matches levels (active=0, passive=1)
% Use numeric context values (match model definition)
context_passive = ones(size(xvals));  % Context = 1 (passive)
context_active  = zeros(size(xvals)); % Context = 0 (active)

% Add dummy AnimalID for random effect group (required)
% dummy_id = repmat("dummy", size(xvals));  % or categorical("dummy")
dummy_id = categorical(repmat("dummy", size(xvals)));

tbl_passive = table(xvals, context_passive, dummy_id, ...
    'VariableNames', {'EngagementProj', 'Context', 'AnimalID'});
tbl_active  = table(xvals, context_active,  dummy_id, ...
    'VariableNames', {'EngagementProj', 'Context', 'AnimalID'});



% Predict using model
[yhat_active,  yci_active]  = predict(lme, tbl_active, 'Conditional', false);
[yhat_passive, yci_passive] = predict(lme, tbl_passive, 'Conditional', false);

%% make plots!

positions = utils.calculateFigurePositions(1, 5, .5, []);

xvals = linspace(min(engagement_proj_all), max(engagement_proj_all), 100);
% Passive line
pred_active = lme.Coefficients.Estimate(1) + ...
               lme.Coefficients.Estimate(2) * xvals;
% Active line
pred_passive = (lme.Coefficients.Estimate(1) + lme.Coefficients.Estimate(3)) + ...
               (lme.Coefficients.Estimate(2) + lme.Coefficients.Estimate(4)) * xvals;
figure(103); clf; hold on;
scatter(engagement_proj_all(context_all==0), sound_proj_all(context_all==0), 5,[0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha',1)
scatter(engagement_proj_all(context_all==1), sound_proj_all(context_all==1), 5,[0.8 0.8 0.8], 'filled', 'MarkerFaceAlpha',1)

% scatter(engagement_proj_all(context_all==0), sound_proj_all(context_all==0), 5,[0.3000    0.2000    0.6000],  'MarkerEdgeAlpha',0.3)
% scatter(engagement_proj_all(context_all==1), sound_proj_all(context_all==1), 5,[0.7000    0.6000    1.000],  'MarkerEdgeAlpha',0.3)
plot(xvals, pred_active, 'k', 'LineWidth', 2.2)
plot(xvals, pred_passive, 'Color',[0.6 0.6 0.6], 'LineWidth', 2.2)


ylabel(ylabel_string)

if nargin > 6
    xlabel(strcat(varargin{1,1},' Projection'))
    ylabel_string_updated = [ylabel_string varargin{1,1}];
else
    xlabel('Engagement Projection')
end

%set figure
set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

% SECOND PLOT WITHOUT ANY SCATTER POINTS
% Plot in figure 103
xvals = linspace(min(engagement_proj_all), max(engagement_proj_all), 100)';

figure(1003);clf; hold on;
fill([xvals; flipud(xvals)], [yci_passive(:,1); flipud(yci_passive(:,2))], ...
     [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4)
fill([xvals; flipud(xvals)], [yci_active(:,1); flipud(yci_active(:,2))], ...
     [0.3 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.3)

plot(xvals, yhat_passive, 'Color',[0.6 0.6 0.6], 'LineWidth', 2)
plot(xvals, yhat_active, 'k', 'LineWidth', 2)
ylabel(ylabel_string)

if nargin > 6
    xlabel(strcat(varargin{1,1},' Projection'))
    ylabel_string_updated = [ylabel_string varargin{1,1}];
else
    ylabel_string_updated = ylabel_string;
    xlabel('Engagement Projection')
end
%set figure
set(gca, 'FontSize', 8, 'Units', 'inches', 'Position', positions(1, :));
utils.set_current_fig;

% Save results
if ~isempty(save_dir)
    mkdir(save_dir)
    cd(save_dir)
    saveas(103,strcat('me_regression_',num2str(ylabel_string_updated),'.fig'));
    exportgraphics(figure(103),strcat('me_regression_',num2str(ylabel_string_updated),'.pdf'), 'ContentType', 'vector');

    saveas(1003,strcat('me_lineonly_regression_',num2str(ylabel_string_updated),'.fig'));
    exportgraphics(figure(1003),strcat('me_lineonly_regression_',num2str(ylabel_string_updated),'.pdf'), 'ContentType', 'vector');

end
