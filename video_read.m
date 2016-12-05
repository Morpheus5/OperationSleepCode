function [video, fps] = video_read(video_file, max_duration)
%VIDEO_READ Loads a video

if ~exist('max_duration', 'var') || isempty(max_duration)
    max_duration = 0;
end

% open video reader
vh = VideoReader(video_file);

% get frame rate
fps = vh.FrameRate;

% store frames
frames = {};

% had frame
while hasFrame(vh)
    % too long
    if max_duration > 0 && vh.CurrentTime > max_duration
        break
    end

    % read frame
    frame = readFrame(vh);

    % add frame
    frames{end + 1} = frame; %#ok<AGROW>
end

% turn into a video
video = cat(1 + ndims(frames{1}), frames{:});

end

