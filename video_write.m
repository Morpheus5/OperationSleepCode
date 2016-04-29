function video_write(video_file, video, fps, format)
%VIDEO_WRITE Write a video file

% infer format from name
if ~exist('format', 'var') || isempty(format)
    switch video_file((end-2):end)
        case 'mp4'
            format = 'MPEG-4';
        case 'm4v'
            format = 'MPEG-4';
        case 'mj2'
            format = 'Motion JPEG 2000';
        case 'avi'
            format = 'Motion JPEG AVI';
        otherwise
            error('Unknown format.');
    end
end

% open writer
vh = VideoWriter(video_file, format);

% set frame rate
if exist('fps', 'var') && ~isempty(fps)
    vh.FrameRate = fps;
end

% set quality
switch format
    case 'Motion JPEG 2000'
        vh.Quality = 80;
    case 'Motion JPEG AVI'
        vh.Quality = 80;
end

open(vh);

% convert format (can not write uint16)
if isa(video, 'uint16') || isa(video, 'single')
    video = im2uint8(video);
end

% write frames
if ndims(video) == 4
    for i = 1:size(video, 4)
        f = im2frame(video(:, :, :, i));
        writeVideo(vh, f);
    end
elseif ndims(video) == 3
    for i = 1:size(video, 3)
        v = repmat(video(:, :, i), 1, 1, 3);
        f = im2frame(v);
        writeVideo(vh, f);
    end
else
    warning('Invalid video passed.');
end

% close
close(vh);

end
