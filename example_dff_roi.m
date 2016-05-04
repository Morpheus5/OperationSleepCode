% loading
load('LNY20RB_RH_9-7-15_SanneSleep6(7).mat');
load('roi/Spontaneous_Data.mat');
load('image_roi/roi_data_image.mat');

% make dff
video_dff = NP_Dff(video.frames(30:end-30));

% draw rois
video_rgb = NP_DrawROIs(video_dff, ROI, squeeze(Spike_Data(1, :, :)));

% save
video_write('dff_LNY20RB_RH_9-7-15_SanneSleep6(7).avi', video_rgb, 30);
