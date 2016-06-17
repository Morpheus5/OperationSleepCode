% clear all; close all;

% Make folders:
mat_dir = 'SpikeVideos'; 
if exist(mat_dir,'dir'); 
    rmdir(mat_dir,'s'); % remove folder with subfolders
end
mkdir(mat_dir);

mat_dir2 = 'DffMoCoVideos'; 
if exist(mat_dir2,'dir'); 
    rmdir(mat_dir2,'s'); % remove folder with subfolders
end
mkdir(mat_dir2);

DIR=pwd; 

% load('image_roi/roi_data_image.mat');

% Make list of videos
mov_listing = dir(fullfile(DIR,'MoCo_MOVIES','*.mat'));
mov_listing = {mov_listing(:).name} % reports which movies it found
data_listing = dir(fullfile(DIR,'MoCoBSBagels_5frameswindow','*.mat'));
data_listing = {data_listing(:).name};

[nblanks formatstring]=progressbar(length(data_listing)); % start progressbar
fprintf(1,['Progress:  ' blanks(nblanks) ]);

% Loop through Dff  (all videos; 1st half of data_listing)
for videoIter = 1:length(mov_listing);
    
    fprintf(1,formatstring,videoIter); % update progressbar
    fprintf(['/' num2str(length(mov_listing)), '\n\n'] );
    tic
    
    [~,file,~] = fileparts(mov_listing{videoIter});
  	load(fullfile(DIR,'MoCo_MOVIES',mov_listing{videoIter}),'video_reg');
    load(fullfile(DIR,'MoCoBSBagels_5frameswindow',data_listing{videoIter}));
    
    savefilename1 = ['SpikeVideos/SpikeROIvideo_' data_listing{videoIter} '.avi'];
    savefilename2 = ['SpikeVideos/SpikeROIvideo_' data_listing{videoIter}];
    savefilename3 = ['DffMoCoVideos/DffBS_' mov_listing{videoIter} '.avi'];
    savefilename4 = ['DffMoCoVideos/DffBS_' mov_listing{videoIter}];
    
    % make dff video:
    video_dff = NP_Dff(video_reg(:,:,30:end-30), 'clim', [0 1]);
    % save video (avi & mat):
    video_write(savefilename3, video_dff, 30);
    video_dff = im2uint16(video_dff);
    save(savefilename4,'video_dff','-v7.3');
    
    toc
    fprintf(1,['\nDffVideo\t', num2str(videoIter), '/', num2str(length(mov_listing)), '\n\n'] );
        
    % draw circles around active ROIs based on DffBS spikes:
    video_rgb = NP_DrawROIs(video_dff, ROI, squeeze(Spike_Data(1, :, :)));
    % save video with spikes (avi & mat):
    video_write(savefilename1, video_rgb, 30);
    video_rgb = im2uint16(video_rgb);
    save(savefilename2,'video_rgb','-v7.3');
    
    toc
    fprintf(1,['\nSpikeVideoDff\t', num2str(videoIter), '/', num2str(length(mov_listing)), '\n\n'] );
    
    clear('Spike_Data')
    load(fullfile(DIR,'MoCoBSBagels_5frameswindow',data_listing{(videoIter+length(mov_listing))}));
    savefilename1 = ['SpikeVideos/SpikeROIvideo_' data_listing{(videoIter+length(mov_listing))} '.avi'];
    savefilename2 = ['SpikeVideos/SpikeROIvideo_' data_listing{(videoIter+length(mov_listing))}];
    
    % draw circles around active ROIs:
    video_rgb = NP_DrawROIs(video_dff, ROI, squeeze(Spike_Data(1, :, :)));
    % save video with spikes (avi & mat):
    video_write(savefilename1, video_rgb, 30);
    video_rgb = im2uint16(video_rgb);
    save(savefilename2,'video_rgb','-v7.3');
    
    toc
    fprintf(1,['\nSpikeVideoAnnuli\t', num2str(videoIter), '/', num2str(length(mov_listing)), '\n\n'] );
    
end

fprintf(1,'\n'); % end progressbar
