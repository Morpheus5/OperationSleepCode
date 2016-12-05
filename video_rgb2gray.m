function video_gs = video_rgb2gray(video_rgb)
%VIDEO_RGB2GRAY Convert RGB video to grayscale video

if 4 ~= ndims(video_rgb) || 3 ~= size(video_rgb, 3)
    error('Unexpected video dimensions. Expecting four dimensions with three color channels.');
end

video_gs = zeros(size(video_rgb, 1), size(video_rgb, 2), size(video_rgb, 4), 'like', video_rgb);
for t = 1:size(video_rgb, 4)
    video_gs(:, :, t) = rgb2gray(video_rgb(:, :, :, t));
end

end

