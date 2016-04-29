% loading
load('SanneSleep23_LNY19RB_12-6-15_II_spont(20).mat');
load('roi/Spontaneous_Data.mat');
load('image_roi/roi_data_image.mat');

% make dff
video_dff = NP_Dff(video.frames(30:end-30));

% draw rois
video_rgb = NP_DrawROIs(video_dff, ROI, squeeze(Spike_Data(1, :, :)));

% save
video_write('dff_SanneSleep23_LNY19RB_12-6-15_II_spont(20).avi', video_rgb, 30);
