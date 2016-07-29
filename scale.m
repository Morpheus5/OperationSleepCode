% Run after your image is already plotted in a figure:

% properties
bar_position = [10 10]; % top left corner of scale bar
bar_scale = 0.5; % microns/pixel
bar_length = 100; % 100 microns
bar_color = [0 1 0];
bar_width = 3;

% calculate coordinates
x1 = bar_position(1);
y1 = bar_position(2);
x2 = round(x1 + bar_length / bar_scale);
y2 = y1;

hold on;

% draw line
plot([x1 x2], [y1 y2], 'Color', bar_color, 'LineWidth', bar_width);

% write text
text((x1 + x2) / 2, y1 + 3 * bar_width, sprintf('%d\\mum', bar_length), 'Color', bar_color, 'FontSize', 20, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');

hold off;