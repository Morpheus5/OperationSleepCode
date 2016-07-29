% load video
VidName = 'SanneSleep31_LNY19RB_1minSpont_1218_(1).mat'
load(VidName);

% concatenate into simpler format
video = cat(3, video.frames(30:end-30).cdata); % if unclipped video
%video = cat(3, video1(30:end-30).cdata); % if clipped video (adjust video#)

% calculate plumes
indices = NP_ExtractPlumes(video);

% calculate dff
video_dff = NP_Dff(video);

durations = [];
intervals = [];
    
% for each plume...
for i = 1:size(indices, 1)
    
    % get start and stop indices
    strt = indices(i, 1);
    stop = indices(i, 2);
    
    durations(end+1) = (stop-strt);
    if i ~= size(indices,1);
        intervals(end+1) = indices(i+1,1)-stop;
    end
        
    % save video
    video_write(sprintf('plume%d.avi', i), video_dff(:, :, strt:stop));
    
    % save image
    im = NP_PixelMass(video(:, :, strt:stop));
    imwrite(im, sprintf('plume%d.jpg', i));
end

PlumeIntervalFiles{end+1} = VidName;
PlumeIntervals(end+1) = ((size(video,3)/30)/size(indices,1));
PlumeAllDurations = [PlumeAllDurations durations];
PlumeAllIntervals = [PlumeAllIntervals intervals];