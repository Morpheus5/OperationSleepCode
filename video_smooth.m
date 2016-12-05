function video_smooth = video_smooth(video, time_steps)
%VIDEO_SMOOTH Smooths a video over time
%   Smooths neighboring frames in a video based on parameter time steps.
%   Each returned frame will be a running average of time_step frames in
%   the input video. Supports both grayscale and multichannel videos.
%   
%   This function takes advantage of matlab n-dimensional convolution
%   function to use a simple smoothing filter along the time dimension of
%   the video. To avoid zero padding, the first and last frames are
%   repeated.

% make filter
filter = ones(1, time_steps, 'like', video) / cast(time_steps, 'like', video);

% amount of padding
if 1 == mod(time_steps, 2)
    pad_before = (time_steps - 1) / 2;
    pad_after = pad_before;
else
    pad_before = (time_steps / 2);
    pad_after = pad_before;
end

% handle different video types
if ndims(video) == 4
    % reshape to filter over time
    filter = reshape(filter, 1, 1, 1, []);
    
    % replicate first and last frames
    video = cat(4, ...
        repmat(video(:, :, :, 1), 1, 1, 1, pad_before), ...
        video, ...
        repmat(video(:, :, :, end), 1, 1, 1, pad_after));
    
    % smooth video
    video_smooth = convn(video, filter, 'valid');
elseif ndims(video) == 3
    % reshape to filter over time
    filter = reshape(filter, 1, 1, []);
    
    % replicate first and last frames
    video = cat(3, ...
        repmat(video(:, :, 1), 1, 1, pad_before), ...
        video, ...
        repmat(video(:, :, end), 1, 1, pad_after));
    
    % smooth video
    video_smooth = convn(video, filter, 'valid');
else
    error('Invalid number of dimensions (expecting 3 or 4).');
end

end

