% load ROI
% load('image_roi/roi_data_image.mat');

% run
% NP_ActivityAroundSpikes(ROI, 'folder_with_raw_video_mat_files', 'folder_to_place_output');
NP_ActivityAroundSpikes(ROI, 'mat', '2secActAroundSpikes', 'window', [1 1]);

