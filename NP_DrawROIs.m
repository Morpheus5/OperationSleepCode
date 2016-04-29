function video_rgb = NP_DrawROIs(frames, roi, spike_data, Fs)
%DRAWROIS Annotate video with regions of interest at spikes
%   For each spike in the spike data, the video will be annotated with a
%   circe indicating the region of interest. This allows for visual
%   inspection of the video to ensure that spike detection matches
%   expectations.

% default parameters
if ~exist('Fs', 'var') || isempty(Fs)
    Fs = 30;
end

% parameters used when drawing
draw_for = round(0.25 * Fs); % 0.25 seconds

% examine inputs
video_length = length(frames);
roi_count = length(roi.stats);

% check lengths
if size(spike_data, 1) ~= roi_count
    error('Expecting %d regions in the spike data.', roi_count);
end
if size(spike_data, 2) ~= video_length
    error('Expecting %d frames in the spike data.', video_length);
end

% generate output video
% concatenate all frames and convert to double)
if isstruct(frames)
    video_gs = cat(3, frames(:).cdata);
else
    video_gs = frames;
end

% convert to RGB
video_r = reshape(video_gs, [], video_length);
video_g = video_r;
video_b = video_r;

% get dimensions
video_height = size(video_gs, 1);
video_width = size(video_gs, 2);

% get spikes to draw
[spike_roi, spike_tm] = find(spike_data);
spike_count = length(spike_roi);

% helper variables
% empty frame...
[x, y] = meshgrid(1:video_width, 1:video_height);
% for calculating time indices
time_idx = 1:video_length;

% colors...
colors = lines(roi_count);
if isa(video_gs, 'uint8')
    colors = uint8(colors .* (256 - 1));
elseif isa(video_gs, 'uint16')
    colors = uint16(colors .* (256 * 256 - 1));
elseif isa(video_gs, 'single')
    colors = single(colors);
end

clear video_gs;

% for each spike...
for i = 1:spike_count
    cur_roi = spike_roi(i);
    
    % stats
    center = roi.stats(cur_roi).Centroid;
    radius2 = (roi.stats(cur_roi).Diameter / 2) ^ 2;
    
    % calculate distance
    distance = (x - center(1)) .^ 2 + (y - center(2)) .^ 2;
    mask_space = distance > radius2 & distance < radius2 * 1.5;
    
    % concatenate mask over time
    mask_time = time_idx >= spike_tm(i) & time_idx <= spike_tm(i) + draw_for;
    
    % set color
    video_r(mask_space, mask_time) = colors(cur_roi, 1);
    video_g(mask_space, mask_time) = colors(cur_roi, 2);
    video_b(mask_space, mask_time) = colors(cur_roi, 3);
end

% combine channels
video_rgb = cat(3, ...
    reshape(video_r, video_height, video_width, 1, video_length), ...
    reshape(video_g, video_height, video_width, 1, video_length), ...
    reshape(video_b, video_height, video_width, 1, video_length));

end

