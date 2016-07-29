function plot_rastergram(tm, spikes, varargin)
%PLOT_RASTERGRAM Plots a rastergram composed of horizontal bars at spikes
%   Requires two parameters:
%
%   tm: Either a vector of times corresponding with the number of columns
%       in spikes, or a scalar representing the same rate of the data
%       (frequency).
%
%   spikes: A logical matrix with rows corresponding with rows in the plot
%           and columns representing time. Each true will be plotted as a
%           single vertical line in the rastergram.

%% parameters

line_width = 1;
colors = [0 0 0];
labels = 1:size(spikes, 1);

% load custom parameters
nparams = length(varargin);
if 0 < mod(nparams, 2)
    error('Parameters must be specified as parameter/value pairs');
end
for i = 1:2:nparams
    nm = lower(varargin{i});
    if ~exist(nm, 'var')
        error('Invalid parameter: %s.', nm);
    end
    eval([nm ' = varargin{i+1};']);
end

%% setup

num_rw = size(spikes, 1); % rows / ROIs
num_tm = size(spikes, 2); % times
num_clrs = size(colors, 1);

% accept frequency
if isscalar(tm)
    tm = (0:(num_tm - 1)) ./ tm;
end

% accept numeric labels
if isnumeric(labels)
    old = labels;
    labels = cell(1, length(old));
    for i = 1:length(old)
        labels{i} = num2str(old(i));
    end
end

% checks
if ~islogical(spikes)
    warning('Spikes should be a logical matrix. Converting...');
    spikes = logical(spikes);
end
if length(tm) ~= num_tm
    error('Vector tm and the number of columns in the spike data must match.');
end
if length(labels) ~= num_rw
    error('The number of labels and the number of rows in the spike data must match.');
end
if size(colors, 2) ~= 3
    error('The color data must have three columns.');
end

% setup plot
xlim([tm(1) tm(end)]);
ylim([0.5 num_rw + 0.5]);
xlabel('Time (s)');
if isempty(labels)
    % no labels
    set(gca, 'YTick', []);
else
    set(gca, 'YTick', 1:num_rw, 'YTickLabel', labels(end:-1:1));
end

% do actual plotting
hold on;
for i = 1:num_rw
    color = 1 + mod(i - 1, num_clrs);
    ln_x = bsxfun(@times, tm(spikes(i, :)), [1; 1; nan]);
    ln_y = repmat(1 + num_rw - i + [0.5; -0.5; nan], 1, size(ln_x, 2));
    plot(ln_x(:), ln_y(:), 'LineWidth', line_width, 'Color', colors(color, :));
end
hold off;

end

