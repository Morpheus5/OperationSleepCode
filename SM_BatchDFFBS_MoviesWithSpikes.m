function SM_BatchDFFBS_MoviesWithSpikes(DIR)

% WALIII
% 09.05.15
% adjusted for Sanne 10.06.15

% This script first downsamples videos, then subtracts the background  
% through disk filter blurring, then blows up video to normal size again  
% for ROI selection later.
% It saves background subtracted videos in mat & avi format.
% It also motion-corrects the BS videos and the data matrix, and makes  
% maximum projection images for ROI selection later.
% Motion correction not working

%% VARIABLES:
mat_dir = 'NewDFFBS_Movies_WithSpikes';

%% Make DFF_MOVIES folder
if exist(mat_dir,'dir'); 
    rmdir(mat_dir,'s'); % remove folder with subfolders
end
mkdir(mat_dir);
if nargin<1 | isempty(DIR); 
    DIR=pwd; 
end

%% Read all videos (in .mat format) and ROI and Spike_Data matrices
mov_listing = dir(fullfile(DIR,'*.mat'));
mov_listing = {mov_listing(:).name} % reports which movies it found

load('roi/Spontaneous_Data.mat');
load('image_roi/roi_data_image.mat');


%% Loop through videos and make dff videos with ROIs indicating spikes
for videoIter = 1:length(mov_listing);
    tic
    [~,file,~] = fileparts(mov_listing{videoIter});
  	load(fullfile(DIR,mov_listing{videoIter}),'video');

    % make dff
    video_dff = NP_Dff(video.frames(30:end-30));

    % draw rois
    video_rgb = NP_DrawROIs(video_dff, ROI, squeeze(Spike_Data(videoIter, :, :)));

    % save
    save_filename = [ fullfile(mat_dir,file) '.avi' ];
    video_write(save_filename, video_rgb, 30);

end

end

