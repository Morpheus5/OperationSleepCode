% load ROI
load('image_roi/roi_data_image.mat');

% run
% NP_ActivityAroundSpikes(ROI, 'folder_with_raw_video_mat_files', 'folder_to_place_output');
NP_ActivityAroundSpikes(ROI, 'clippedmovies', '2secActAroundSpikesNew2', 'window', [1 1]);

% To generate the PixelMass figure:
% 
% figure;
% imshow(img_comp);
% colormap(cmap); colorbar;