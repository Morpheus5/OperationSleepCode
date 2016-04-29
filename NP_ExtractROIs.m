function roi_ave = NP_ExtractROIs(frames, roi, donut)
%NP_EXTRACTROIS Extract regions of interest
%   This function takes a video (`frames`) and a bunch of regions of
%   interest (`roi`), and will extract and average trace for each ROI.
%
%   Optionally, a `donut` multiplier can be specificied. If the donut is
%   above zero, a donut (annulus) is created around the region of interest
%   and subtracted from the region of interest time course. The variable
%   donut specifices how many pixels should be in the donut. For example,
%   if `donut` = 1, then the annulus will contain the same number of pixels
%   as the region of interest.

% default value is no donut
% use donut = 2 for roughtly equal number of pixels in donut and donut hole
if ~exist('donut', 'var') || isempty(donut)
    donut = 0;
end

% examine inputs
video_length = length(frames);
roi_count = length(roi.stats);

% convert to grayscale video
if isstruct(frames)
    video_gs = cat(3, frames(:).cdata);
else
    video_gs = frames;
end

% get roi data

% make grid (number each pixel of the video)
[x, y] = meshgrid(1:size(video_gs, 2), 1:size(video_gs, 1));

% allocate memory for return matrix (1 row per ROI, 1 column per frame)
roi_ave = zeros(roi_count, video_length);

% reshape video to have one row per pixel
video_gs_r = reshape(video_gs, [], video_length);

% for each region of interest...
for i = 1:roi_count
    % get center and radius squared
    center = roi.stats(i).Centroid;
    radius2 = (roi.stats(i).Diameter / 2) ^ 2;
    
    % calculate distances between pixels and ROI center
    d = (x - center(1)) .^ 2 + (y - center(2)) .^ 2;
    
    if 0 < donut
        % has donut, should subtract it
        % make masks representing donut and donut hole
        mask_donut_hole = d < radius2;
        mask_donut = d >= radius2 & d < (radius2 * (1 + donut));
        
        % get donut and donut hole
        ave_donut_hole = mean(video_gs_r(mask_donut_hole(:), :), 1);
        ave_donut = mean(video_gs_r(mask_donut(:), :), 1);
        
        % subtract donut from donut hole
        roi_ave(i, :) = ave_donut_hole - ave_donut;
    else
        % make mask
        mask = d < radius2;
    
        % extract mean time series
        roi_ave(i, :) = mean(video_gs_r(mask(:), :), 1);
    end
end

end

