function video_bs = NP_Bs(frames)

% WALIII
% 09.05.15
% adjusted for Sanne 10.06.15

% Nathan
% made more generic, separated from batching functionality 4/28/16

% This script first downsamples videos, then subtracts the background  
% through disk filter blurring, then blows up video to normal size again  
% for ROI selection later.
% It saves background subtracted videos in mat & avi format.
% It also motion-corrects the BS videos and the data matrix, and makes  
% maximum projection images for ROI selection later.
% Motion correction not working

%% VARIABLES:
downsample_factor = 1; %downsampling: 4; otherwise 1
motion_correction = false;
use_gpu = gpuDeviceCount > 0;

%% Downsample, subtract background and smooth
video_bs = cat(3, frames(:).cdata);

if use_gpu
    video_bs = gpuArray(video_bs);
end

% downsampling: downsample_factor = 0.25; otherwise 1
if downsample_factor ~= 1
    video_bs = imresize(video_bs, 1 / downsample_factor);
end

% disk size for blurring to subtract from data for background subtraction; 
filter_disk_blur = fspecial('disk', 50 / (downsample_factor ^ 2));
bground = imfilter(video_bs, filter_disk_blur, 'replicate');
video_bs = video_bs - bground;

filter_disk_smooth = fspecial('disk', 1); % smoothing with a very small disk
video_bs = imfilter(video_bs, filter_disk_smooth);

%% scale back to size after downsampling
if downsample_factor ~= 1
    video_bs = imresize(video_bs, downsample_factor);
end

if use_gpu
    video_bs = gather(video_bs);
end

%% correct motion artifacts in movies 
% Commented out for now - veeerrrry slow, sometimes even crashes
if motion_correction
    error('Not yet supported. Must test.');

    [optimizer, metric] = imregconfig('multimodal'); % Same device, but might have different brightness ranges

    for frameIter=1:(size(video_bs, 3)); 
        video_bs(:,:,frameIter)=imregister(video_bs(:, :, frameIter), video_bs(:,:,100), 'rigid', optimizer, metric);
    end
end

%% Rescale
video_bs = video_adjust(video_bs, [0.5 0.999]);

end
