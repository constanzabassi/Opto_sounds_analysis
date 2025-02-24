function place_text_labels(labels, colors, y_offset_base, fontSize)
    % Get current axis limits
    x_range = xlim;
    y_range = ylim;
    % Calculate base text position
    text_x = x_range(2) - y_offset_base * diff(x_range);
    text_y = y_range(1) + y_offset_base * diff(y_range);
    % Default font size if not provided
    if nargin < 4
        fontSize = 12;
    end
%     % Ensure labels and colors match
%     if length(labels) ~= size(colors,1)
%         error('Labels and colors must have the same length');
%     end
    % Auto-calculate evenly spaced y-offsets
    num_labels = length(labels);
    y_offsets = linspace(0, 0.2 * (num_labels - 1), num_labels); % Adjusted scaling
    % Place text labels
    for i = 1:num_labels
        text(text_x, text_y - y_offsets(i) * diff(y_range), labels{i}, ...
             'Color', colors(i,:), 'FontSize', fontSize);
    end
end