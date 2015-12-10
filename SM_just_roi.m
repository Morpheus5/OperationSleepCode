function [roi_ave] = SM_just_roi(ROI,varargin)
% ROIS is a cell array of image indices returned by FS_image_roi

% VARIABLES
ave_fs=30; % frame rate
duration=30; % movie duration in seconds
%ave_fs_i=20*ave_fs; % after interpolation
resize = 1;
colors=eval(['winter(' num2str(length(ROI.coordinates)) ')']);
sono_colormap='hot';
save_dir='roi';
template=[];

%%
% first convert ROIS to row and column indices

if resize~=1
	disp(['Adjusting ROIs for resizing by factor ' num2str(resize)]);

	for VideoIter=1:length(ROI.coordinates)
		ROI.coordinates{VideoIter}=round(ROI.coordinates{VideoIter}.*resize);
	end
end

mkdir(save_dir);
mov_listing=dir(fullfile(pwd,'*.mat'));
mov_listing={mov_listing(:).name};

to_del=[];
for VideoIter=1:length(mov_listing)
	if strcmp(mov_listing{VideoIter},'dff_data.mat')
		to_del=VideoIter;
	end
end
mov_listing(to_del)=[];

%%
% Make matrix for video data (mov_data)
% Here, to calculate matrix size; below for data when looping through
% videos
load(fullfile(pwd,mov_listing{1}),'video');
mov_data_temp = video.frames;
for VideoIter = 1:length(mov_data_temp)
    mov_data(:,:,VideoIter) = mov_data_temp(VideoIter).cdata;
end
mov_data = double(mov_data);

[rows,columns,frames]=size(mov_data);

t = 1/ave_fs % seconds/frame
ave_time = (t):(t):(duration);
number_of_frames = length(ave_time)
roi_ave.raw={};

clear mov_data % to allow looping through each video below

%% loop through videos, extract frames
disp('Looping through each video...');
Videos=length(mov_listing);
for VideoIter=1:Videos
	disp(['Processing video ' num2str(VideoIter) ' of ' num2str(length(mov_listing))]);
	load(fullfile(pwd,mov_listing{VideoIter}),'video');
    mov_data_temp = video.frames;
 
    %% loop through frames, extract row and column values
        for FrameIter = 1:length(mov_data_temp)  
             mov_data(:,:,FrameIter) = mov_data_temp(FrameIter).cdata(:,:); % mov_data=x512,y512,frameiter
        end
        mov_data = double(mov_data); 
        
    %%
        [~,file,~]=fileparts(mov_listing{VideoIter}); % [path,file,extension]
        save_file=[ file '_roi' ];

        [~,~,frames] = size(mov_data); %[rows,columns,frames]
        roi_n = length(ROI.coordinates);
        roi_data = zeros(roi_n,frames);    
            
    %% extract mean values per frame in each ROI circle              
        disp('Computing ROI averages...');

        [nblanks formatstring] = fb_progressbar(100);
        fprintf(1,['Progress:  ' blanks(nblanks)]);

        % unfortunately we need to for loop by frames, otherwise
        % we'll eat up too much RAM for large movies

        for ROIiter = 1:roi_n
            fprintf(1,formatstring,round((ROIiter/roi_n)*100));

            for FrameIter2 = 1:frames 
                roi_data_tmp = mov_data(ROI.coordinates{ROIiter}(:,2),ROI.coordinates{ROIiter}(:,1),FrameIter2); 
                roi_data(ROIiter,FrameIter2,:) = mean(roi_data_tmp(:)); 
            end
        end

        fprintf(1,'\n');
        
% %% interpolate mean values in each ROI
%         for ROIiter2=1:roi_n
%             roi_t = roi_data(ROIiter2,:,:);% ROIs,frames,mean values 
%             
%             frame_idx = 1:length(roi_t);
%             timevec=(frame_idx./30);           
%             interp_raw=interp1(timevec,roi_t,ave_int_time,'spline');
%             
%         end
%       
%% Save files        
    %roi_ave.mov_data{VideoIter}=mov_data; % This line is untested, this should save the actual raw data (unaveraged)
    roi_ave.raw{VideoIter}=(roi_data); 
%     roi_ave.interp_raw(ROIiter2,:,VideoIter) = interp_raw;

    roi_ave.filename{VideoIter}=mov_listing{VideoIter};
    roi_ave.t=ave_time;
%     roi_ave.int_=ave_int_time;
    save(fullfile(save_dir,['ave_roi.mat']),'roi_ave');
    disp('Generating roi_ave matrix...');

end