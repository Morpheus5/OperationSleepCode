function FS_BatchDff_NEW(DIR, varargin)

%run thorough directory and make Background Subtracted Movies in AVI format
% This Script is for 'unprocessed' videos


% WALIII
% For unprocessed videos
% 09.05.15

filt_rad=1; % gauss filter radius
filt_alpha=1; % gauss filter alpha
lims=3; % contrast prctile limits (i.e. clipping limits lims 1-lims)
cmap=colormap('jet');
per=10; % baseline percentile (0 for min)
counter = 1;

mat_dir='DFF_MOVIES';
counter = 1;

if exist(mat_dir,'dir') rmdir(mat_dir,'s');
end
mkdir(mat_dir);

outlier_flag=0;
if nargin<1 | isempty(DIR), DIR=pwd; 
end
mov_listing=dir(fullfile(DIR,'*.mat'));
mov_listing={mov_listing(:).name};
filenames=mov_listing;

disp('Creating Dff movies');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(mov_listing)
clear video; 
clear NormIm; 
clear dff2; 
clear dff; 
clear mov_data3; 
clear mov_data2; 
clear mov_data;
clear mov_data_16bit;
clear baseline;
clear test;
clear v;

[path,file,ext]=fileparts(filenames{i});
fprintf(1,formatstring,round((i/length(mov_listing))*100));

load(fullfile(DIR,mov_listing{i}),'video')
save_filename = [ fullfile(mat_dir,file) ];
LastFrame = length(video.frames);
mov_data = video.frames;


%mov_data_16bit = zeros(512,512,LastFrame);
for frameIter = 1:LastFrame-3;
    mov_data_16bit1 = uint16(mov_data(frameIter).cdata);
    mov_data_16bit2 = uint16(mov_data(frameIter+1).cdata);
    mov_data_16bit3 = uint16(mov_data(frameIter+2).cdata);
    mov_data_16bit(:,:,frameIter) = uint16(mov_data_16bit1 + mov_data_16bit2 + mov_data_16bit3)/3;
end
clear mov_data_16bit1; clear mov_data_16bit2; clear mov_data_16bit3;


test = single(mov_data_16bit(:,:,12:end));
[rows,columns,frames]=size(test);
  
%%%=============[ FILTER Data ]==============%%%

disp('Gaussian filtering the movie data...');

h=fspecial('gaussian',filt_rad,filt_alpha);
test=imfilter(test,h,'circular');

disp(['Converting to df/f using the ' num2str(per) ' percentile for the baseline...']);

baseline=repmat(prctile(test,per,3),[1 1 frames]);

 h=fspecial('gaussian',20,40);
 baseline = imfilter(baseline,h,'circular'); % filter baseline

tot = (test-baseline); 
%baseline(baseline<0)=1;
dff = (tot./(baseline)).*100;
dff2 = imresize(dff,1);% Scale Data

H = prctile(mean(max(dff2(:,:,15:end))),99);
L = prctile(mean(mean(dff2(:,:,15:end))),50);
    
clim = [double(L) double(H)];
    
NormIm(:,:,:) = mat2gray(dff2, clim);

%figure(1); for  iii = 7:size(NormIm,3);  IM(:,:) = NormIm(:,:,iii); imagesc(IM); pause(0.05); end



%% Write VIDEO

v = VideoWriter(save_filename);
v.Quality = 80;
v.FrameRate = 30;

open(v)

for ii = 2:size(NormIm,3); 
figure(1);
colormap(gray);
IM(:,:) = NormIm(:,:,ii);
writeVideo(v,IM)
imagesc(IM); 
end
close(v)

%% Save Data from aggregate
% Test = TotalX2;
%mov_data = video.frames;
%im_resize = 1;

% save(save_filename,'test','mov_data','im_resize','-v7.3')
end