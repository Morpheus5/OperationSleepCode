% make some data
% rows represent rows in the plot, columns represent time
data = rand(20, 10000) > 0.995;
data([1 2 5 7 8 9 10], 2000) = true;
data(:, 5000) = true;

tm = (1:10000) ./ 250;

plot_rastergram(tm, data, 'colors', lines(20));