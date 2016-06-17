function video = video_register(video, reference, show_progress)
%VIDEO_REGISTER

if 3 ~= ndims(video)
    error('Expecting a grayscale video.');
end

% show progress
if ~exist('show_progress', 'var')
    show_progress = true; % default value
end

% show progress
if show_progress
    h = waitbar(0, 'Registering video...');
end

% convert
if isa(video, 'double')
    convert = false;
else
    convert = class(video);
    video = double(video);
end

% reference frame
if ~exist('reference', 'var') || isempty(reference)
    % DEFAULT: middle frame
    reference = -1;
    ref_frame = video(:, :, round(size(video, 3) / 2));
elseif reference == 0
    % mean
    ref_frame = mean(video, 3);
else
    ref_frame = video(:, :, reference);
end

ref_fft = fft2(ref_frame);

for i = 1:size(video, 3)
    [~, out_fft] = dftregistration(ref_fft, fft2(video(:, :, i)));
    video(:, :, i) = abs(ifft2(out_fft));
    
    % update progress
    if show_progress
        waitbar(i / size(video, 3));
    end
end

% close progress
if show_progress
    close(h);
end

% undo convert
if convert
    video = cast(video, convert);
end

end
