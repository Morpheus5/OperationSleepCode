function [mov_data_BS] = SM_BatchBS(DIR)

% WALIII
% 09.05.15
% adjusted for Sanne 10.06.15

% This script first downsamples videos, then subtracts the background  
% through disk filter blurring, then blows up video to normal size again  
% for ROI selection later.
% It saves background subtracted videos in mat & avi format.
% It also motion-corrects the BS videos and the data matrix, and makes  
% maximum projection images for ROI selection later.

%% VARIABLES:
DSFactor = 4; %downsampling: 4; otherwise 1
mat_dir = 'BS_MOVIES';

%% Make DFF_MOVIES folder
if exist(mat_dir,'dir'); 
    rmdir(mat_dir,'s'); % remove folder with subfolders
end
mkdir(mat_dir);
if nargin<1 | isempty(DIR); 
    DIR=pwd; 
end

%% Read all videos (in .mat format)
mov_listing = dir(fullfile(DIR,'*.mat'));
mov_listing = {mov_listing(:).name} % reports which movies it found

%% Preallocate matrices
mov_data_16bit = zeros(512,512,length(mov_listing));
TotalMaX = zeros(512,512,length(mov_listing));

%% Loop through videos
for videoIter = 1:length(mov_listing);
    [~,file,~] = fileparts(mov_listing{videoIter});
  	load(fullfile(DIR,mov_listing{videoIter}),'video');
    MaxImages = zeros(512,512,length(video.frames)); % preallocate matrix 
    mov_data_BS_MC = zeros(512,512,length(video.frames)); % preallocate matrix 
        % BS = background subtracted; MC = motion corrected
    
        % open window for video
        figure, set(gcf, 'Color','white');
        axis tight;
        set(gca, 'nextplot','replacechildren', 'Visible','off');

        % create AVI object
        save_filename = [ fullfile(mat_dir,file) ];
        vidObj = VideoWriter(save_filename);
        vidObj.Quality = 30;
        vidObj.FrameRate = 30;
        open(vidObj);
        colormap(bone);
                       
        % Read in video and save in new matrix
        mov_data = video.frames;
        clear video.frames; % so that after this loop, the next video can 
        % be read;   mov.data matrix = videoframes X, Y, frames
              
        %% Downsample, subtract background and smooth  
        for frameIter = 1:(length(mov_data));
           mov_data_16bit(:,:,frameIter) = uint16(mov_data(frameIter).cdata/256);
        end

        mov_data_DS = imresize((mov_data_16bit),(1/DSFactor)); 
        % downsampling: DSFactor = 0.25; otherwise 1

        DiskBlur = fspecial('disk',(80/(DSFactor^2))); 
        % disk size for blurring to subtract from data for background subtraction; 
        bground = imfilter(mov_data_DS,DiskBlur);
        mov_data_BS = mov_data_DS - bground;

        DiskSmooth = fspecial('disk',1); % smoothing with a very small disk
        mov_data_BS = imfilter(mov_data_BS,DiskSmooth);

        %% Calculating reference frame for imregister (100th frame) 
        LinKatOne = 1;
        for u = 50:55; % take y-lines in the middle of the downsized image
        LinKatTemp = cat(1,mov_data_BS(:,u,100)); % for 100th frame
        LinKat1 = cat(1,LinKatOne,LinKatTemp); 
        end
        % makes one large vector to later go through the Fl values for scaling
        for Y_iter = 2:size(mov_data_BS,2);
        LinKatTemp2 = cat(1,mov_data_BS(:,Y_iter,100));
        LinKat2 = cat(1,LinKat1,LinKatTemp2);
        end
        
        H = prctile(LinKat2,95)+10; % determine percentile values for
        L = prctile(LinKat2,5);      % image & save so I can check later
         
        save('LK2','LinKat2')
        
        %% scale back to size after downsampling
        mov_data_BS=imresize(mov_data_BS,DSFactor);
                            
        %% correct motion artifacts in movies 
        % Commented out for now - veeerrrry slow, sometimes even crashes
                        
%         [optimizer, metric] = imregconfig('multimodal'); % Same device, but might have different brightness ranges
%         
%         startMotionCorr = toc
%         
%         for frameIter=1:(size(mov_data_BS,3)); 
%             mov_data_BS_MC(:,:,frameIter)=imregister(mov_data_BS(:,:,frameIter),mov_data_BS(:,:,100),'rigid',optimizer,metric);
%         
%             figure(1);
%             colormap(bone);
%             image(mov_data_BS_MC(:,:,frameIter),'CDataMapping','scaled');
%             caxis([double(L),double(H)]); 
%             set(gca,'ydir','reverse'); % Otherwise the y-axis would be flipped
%             writeVideo(vidObj, getframe(gca));
%         end
        close(gcf); % gcf = get current figure 
        close(vidObj);
    
        %% make maximum projection images of each video
        
        MaxProj = max(mov_data_BS_MC,[],3);
        colormap(bone);
        imagesc(MaxProj);
        maX = mat2gray(MaxProj);
        maX = im2uint16(maX);
        imwrite(maX,save_filename,'png');

        TotalMaX(:,:,videoIter) = maX;
       
        %% Save background subtracted videos in BS_MOVIES folder
        for frameIter = 1:size(mov_data_BS_MC,3);
            video.frames(frameIter).cdata(:,:,:) = mov_data_BS_MC(:,:,frameIter);
            video.frames(frameIter).colormap = [];
            % matrix design: video.frames{frames}.cdata(X,Y,uint16)
        end
        save(save_filename,'video','-v7.3');

end

%% create maximum projection images with corrected alignment 
[optimizer, metric] = imregconfig('multimodal');
for frameIter = 1:size(TotalMaX,3);
    MaxImages(:,:,frameIter) = imregister(TotalMaX(:,:,frameIter), TotalMaX(:,:,1),'rigid',optimizer, metric);
end;
 
%% create a maximum image of all videos to use for ROI selection later 
Max_MaxProj = max(MaxImages,[],3);  % Was totalX, but I think you'd want  
                                    % to take tiledImage

imwrite(Max_MaxProj,'Max_BSimage','png');

end
