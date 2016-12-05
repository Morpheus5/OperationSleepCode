% load comparison data
load('../Comparison/ROI_data_cleansed.mat');
load('../Comparison/roi_spatial_stats.mat');

% read in video
video_gs = video_read('../Comparison/2016-01-21 07 06 34.mov');
video_gs = video_rgb2gray(video_gs);

% convert ROI format
% roi_stats(5) corresponds with 2016-01-21
ROI = struct('stats', roi_stats(5).bw_stats);

% extract roi_ave
donut = 1;
roi_trace = NP_ExtractROIs(video_gs(:, :, 5:end), ROI, donut);

% show traces
figure;
n = size(roi_trace, 1); % number of regions
t = (1:size(roi_trace, 2)) ./ 30; %times
clrs = lines(n); % make colorful lines
subplot(1, 3, 1); % plot image for reference
im_summary = mean(video_gs, 3);
imshow(imadjust(im_summary ./ max(im_summary(:))));
for i = 1:n
    h = viscircles(ROI.stats(i).Centroid, ROI.stats(i).Diameter / 2, 'Color', clrs(i, :));
end
subplot(1, 3, 2); % plot area for traces
plot_many(t, roi_trace', clrs);
subplot(1, 3, 3); % plot spikes
[~, ~, spikes] = NP_From_Raw_To_Traces({roi_trace}, 2);
%[~, ~, spikes] = Byrons_From_Raw_To_Traces({roi_trace});
plot_many(squeeze(spikes(1, :, :))');

% make dff
video_dff = NP_Dff(video_gs(:, :, 5:end));

% annotate ROIS
video_rgb = NP_DrawROIs(video_dff, ROI, squeeze(spikes(1, :, :)));
if 0 < donut
    video_write('donut_dff_2016-01-21 07 06 34.avi', video_rgb, 30);
    print(gcf, 'donut_trace_2016-01-21 07 06 34.avi', '-dpng', '-r300');
else
    video_write('dff_2016-01-21 07 06 34.avi', video_rgb, 30);
    print(gcf, 'trace_2016-01-21 07 06 34.avi', '-dpng', '-r300');
end
