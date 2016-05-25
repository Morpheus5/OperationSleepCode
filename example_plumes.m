% load video
load('../Recordings/SanneSleep23_LNY19RB_12-6-15_II_spont(20).mat');

% concatenate into simpler format
video = cat(3, video.frames(30:end-30).cdata);

% calcualte plumes
indices = NP_ExtractPlumes(video);

% calculate dff
video_dff = NP_Dff(video);

% for each plume...
for i = 1:size(indices, 1)
    % get start and stop indices
    strt = indices(i, 1);
    stop = indices(i, 2);
    
    % save video
    video_write(sprintf('plume%d.avi', i), video_dff(:, :, strt:stop));
    
    % save image
    im = NP_PixelMass(video(:, :, strt:stop));
    imwrite(im, sprintf('plume%d.jpg', i));
end
