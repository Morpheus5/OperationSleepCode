% load video
% load Spontaneous_Data

% concatenate into simpler format
video = cat(3, video.frames(30:end-30).cdata); % if unclipped video
%video = cat(3, video1(30:end-30).cdata); % if clipped video (adjust video#)


videoname = 'Video(1)';
window_before = 30; 
window_after = 30;
total_length = window_before + 1 + window_after; 
fs = 30;

% Maybe add something here were I can select a certain cell

% New matrices -- select ROI of interest
cell18 = find(Spike_Data(1,18,:));
extraction18 = zeros(length(cell18), 2);
roi_trace = squeeze(Trace_Data(1,:,:));

% Save all spiking indices: 
indices = zeros(length(cell18),2);
for i = 1:length(cell18)
    indices(i,1) = (cell18(i) - window_before);
    indices(i,2) = (cell18(i) - window_after);
end

num_spikes = size(indices, 1);
num_frames = indices(1, 2) - indices(1, 1);

% prepare matrices to save all traces
sum_traces = zeros(num_spikes, num_frames);
images = [];

% for each spike...
for i = 1:num_spikes
    
    % get start and stop indices
    strt = indices(i, 1);
    stop = indices(i, 2);
    
    % Check if windows don't overlap:
    
    % Make sum_traces matrix - here for ROI 18
    sum_traces(i,:) = Trace_Data(1, 18, strt:stop); % adjust for ROI#
    
    % save video
    video_write(sprintf('spike%d.avi', i), video(:, :, strt:stop));
    
    % save image
    im = NP_PixelMass(video(:, :, strt:stop));
    imwrite(im, sprintf('spike%d.jpg', i));
    
    % Save max projection images:
    StdProj = std(images,0,3); % 0 = no weighting, 3 = std over 3rd dimension of matrix
    %  colormap(bone);
    %  figure(); imagesc(StdProj); 
    maX18 = mat2gray(StdProj);
    maX18 = im2uint16(maX18);
    save_filename = ([videoname '_STD_cell18_' num2str(indices(i,1)) '-' num2str(indices(i,2)) '.png']);
    imwrite(maX18,save_filename,'png');
    
end


