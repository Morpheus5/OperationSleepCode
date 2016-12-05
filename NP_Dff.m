function video_dff = NP_Dff(frames, varargin)
%FS_DFF_NEW Generates delta f over f video
%   When passed a video, this function calculates a delta f over f version
%   of the video, which is useful for visualizing Calcium traces. This
%   process is done in a few steps:
%
%   1. The video is smoothed over time (each frame is averaged with up to
%   three neighboring frames)
%
%   2. A small gaussian filter is applied over the video (spatial
%   smoothing).
%
%   3. A baseline image is calcualted using the 3rd percentile value for
%   each pixel over time. The baseline image is gaussian filtered using a
%   much larger filter (spatial smoothing).
%
%   4. The square of the baseline image is subtracted from the square of
%   the pixel intensity (roughly delta f) and this is then divided by the
%   baseline (roughly delta f over f). The square is used based on Will's
%   feedback that it proves more reliable at suppressing low variance
%   pixels.

% parameters
filt_rad = 1; % gauss filter radius
filt_alpha = 1; % gauss filter alpha
per = 10; % baseline percentile (0 for min)
clim = [0.5 0.99]; % color range

% load custom parameters
nparams = length(varargin);
if 0 < mod(nparams, 2)
    error('Parameters must be specified as parameter/value pairs.');
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
    mov = cat(3, frames(:).cdata);
else
    mov = frames;
end

% smooth (use a [1/3 1/3 1/3] convolution along the third dimension)
mov = video_smooth(single(mov), 3);

% Gaussian filter video
h = fspecial('gaussian', filt_rad, filt_alpha);
mov = imfilter(mov, h, 'replicate');

% calculate baseline using percentile and repeat over video
baseline = prctile(mov, per, 3);

% use circular gaussian filter for baseline
h = fspecial('gaussian', 20, 40);
baseline = imfilter(baseline, h, 'replicate'); % filter baseline

% calculate dff
video_dff = bsxfun(@minus, mov .^ 2, baseline .^ 2);
video_dff = bsxfun(@rdivide, video_dff, baseline);

% get high and low percentiles
video_dff = video_adjust(video_dff, clim);

end

