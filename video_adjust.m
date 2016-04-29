function video = video_adjust(video, tol)
%VIDEO_ADJUST Summary of this function goes here
%   Detailed explanation goes here


if ~exist('tol', 'var') || isempty(tol)
    tol = [0.2 1];
end

% handle non image values gracefully
if isfloat(video)
    mn = min(video(:));
    mx = max(video(:));
    if mn < 0 || mx > 1
        video = (video - mn) ./ (mx - mn);
    end
end

old_dim = size(video);
clim = stretchlim(video(:), tol);
video = imadjust(video(:), clim, []);
video = reshape(video, old_dim);

end

