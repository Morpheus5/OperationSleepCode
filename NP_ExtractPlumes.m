function indices = NP_ExtractPlumes(frames, varargin)
%NP_EXTRACTPLUMES Find periods of video that may be plumes
%   This uses mean intensity per frame to identify potential plumes in the
%   video. The process involves a few key steps:
%
%   1. A window of 300 (customizable) frames is used to calculate a 
%   threshold based on the 90th percentile (customizable) intensity value.
%
%   2. A smoothed intensity curve is then compared to this threshold to
%   determine plumes.
%
%   3. The plume period is then extended forward and backward to include
%   any ramp up and ramp down period, based on the derivative of the
%   smoothed intensity curve.
%
%   The function returns indices that correspond with the start and stop of
%   each plume found. Returned value `indices` is a n x 2 matrix, where
%   each row corresponds to a single plume.

% parameters
moving_avg_window = 300; % frames over which threshold is caculated
moving_avg_ptile = 90; % percentile used for threshold
min_duration = 5; % in frames

% load custom parameters
nparams=length(varargin);
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

% turn into movie data
if isstruct(frames)
    video = cat(3, frames(:).cdata);
else
    video = frames;
end

% calculate mean intensity
intensity = mean(reshape(video, [], size(video, 3)));

% smooth intensity
smooth_intensity = sgolayfilt(intensity, 3, 41);

% calculate threshold
col = (1:moving_avg_window)';
row = 0:(size(video, 3) - moving_avg_window);
idx = bsxfun(@plus, col, row);
threshold = prctile(intensity(idx), moving_avg_ptile);
threshold = [threshold(1) * ones(1, ceil((size(video, 3) - length(threshold)) / 2)) threshold threshold(end) * ones(1, floor((size(video, 3) - length(threshold)) / 2))];

% % visualize
% t = 1:length(intensity);
% figure; plot(t, intensity, t, smooth_intensity, t, threshold);

% find when smooth intensity is above threshold
above_threshold = (smooth_intensity > threshold);

% expand periods when above threshold based on derivative
dx = [0 diff(smooth_intensity)];

% extend backwards
% as long as smoothed intensity is increasing
while true
    should_include = (dx(1:(end-1)) > 0 & ~above_threshold(1:(end-1)) & above_threshold(2:end));
    if any(should_include)
        above_threshold([should_include false]) = true;
    else
        break;
    end
end

% extend forwards
% as long as smoothed intensity is decreasing
while true
    should_include = (dx(2:end) < 0 & ~above_threshold(2:end) & above_threshold(1:(end-1)));
    if any(should_include)
        above_threshold([false should_include]) = true;
    else
        break;
    end
end

% % plot periods above threshold
% t = 1:length(intensity);
% figure; scatter(t(above_threshold), intensity(above_threshold));

% based on:
% http://stackoverflow.com/questions/3274043/finding-islands-of-zeros-in-a-sequence

% calculate differences
dsig = diff([0; above_threshold(:); 0]);
index_start = find(dsig > 0);
index_end = find(dsig < 0) - 1;
duration = index_end - index_start + 1;

% threshold duration
index_string = (duration >= min_duration);
index_start = index_start(index_string);
index_end = index_end(index_string);

% combine indices
indices = [index_start index_end];

end
