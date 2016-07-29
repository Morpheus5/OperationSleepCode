function [] = SM_BatchMoCo(DIR)

mat_dir = 'MoCo_MOVIES'; % MoCo = MotionCorrection

%% Make MoCo_MOVIES folder
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

%% Loop through videos
for videoIter = 1:length(mov_listing);
    [~,file,~] = fileparts(mov_listing{videoIter});
  	v = load(fullfile(DIR,mov_listing{videoIter}),'video*');
    fn = fieldnames(v);
    video = v.(fn{1});
    clear v;
    
    save_filename = [ fullfile(mat_dir,file) ];
    save_filename2 = [ fullfile(mat_dir,file) '.avi'];
    mov_data = cat(3,video(:).cdata);
    mov_data = mat2gray(mov_data);
    
    %% Run DFT-based subpixel alignment for motion correction:
    
    % https://www.osapublishing.org/ol/abstract.cfm?uri=ol-33-2-156
    % Efficient subpixel image registration algorithms
    % Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup
   
    video_reg = video_register(mov_data);
    
    % fix any out of bound values due to video registration
    video_reg(video_reg > 1) = 1;
    video_reg(video_reg < 0) = 0;
    
    % write video:
    video_write(save_filename2, video_reg, 30);
    
    % save mat file:
    video_reg = im2uint16(video_reg);
    save(save_filename,'video_reg','-v7.3');
end
