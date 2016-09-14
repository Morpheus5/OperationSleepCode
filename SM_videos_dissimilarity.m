% load in the videos with activity around spikes  
% saved from NP_ActivityAroundSpikes 

DIR = pwd;
videolist = dir(fullfile(DIR,'*raw.mp4')); 
videolist = {videolist(:).name};
NrSpikes = length(videolist);
n = [1]; % ROIs of interest

% prepare matrices for video dissimilarity scores:
total_array = [];
d = zeros(NrSpikes, NrSpikes);
difference = [];
sumd = [];

% Load in video:
for iter = 1:NrSpikes
    [~,file,~] = fileparts(videolist{iter});
  	video = video_read(fullfile(DIR,videolist{iter}));
    if 4 == ndims(video)
        video = video_rgb2gray(video);
    end
    
    % now that we know the dimensions, setup our matrices
    if iter == 1
        total_array = zeros(NrSpikes, size(video, 1), size(video, 2), size(video, 3));
    end
    
    % add image to array
    total_array(iter, :, :, :) = video;
    
end

% make difference matrix 'd' & summed dissimilarity vector 'sumd':
for i = 1:NrSpikes
   for j = (i+1):NrSpikes
       % subtract videos from each other frame by frame
       temp = total_array(i, :, :, :) - total_array(j, :, :, :);
       
       % save subtracted images:
       % std image:
       diff = std(squeeze(temp),1,3);
       figure(); imagesc(diff); colormap('gray(2555');
       savefilename = (['Subtracted_ROIs' num2str(i) '-' num2str(j) 'Std.jpg']);
       saveas(gcf, savefilename); close(gcf);
       % mean image:
       diff = mean(squeeze(temp),3);
       figure(); imagesc(diff); colormap('gray(2555');
       savefilename = (['Subtracted_ROIs' num2str(i) '-' num2str(j) 'Mean.jpg']);
       saveas(gcf, savefilename); close(gcf);
       
       % Make sumd vector:
       d(i, j) = mean(abs(temp(:))); % MAE: mean absolute error
       sumd(end+1) = d(i, j);
    end
end

% save summed dissimilarity vector for all sleep cells:
save('Dissimilarity_matrix','sumd','-v7.3')

% plot summed dissimilarity distribution:
figure(1); hold on; title('Dissimilarity histogram videos');
hist(sumd, 200); hold off
saveas(1,'Dissimilarity Histogram videos.eps')

%% ROI based approach

DIR = pwd;
traceslist = dir(fullfile(DIR,'2016*.mat')); 
traceslist = {traceslist(:).name};
NrSpikes = length(traceslist);

% prepare matrices for video dissimilarity scores:
traces_sumd = [];

% Load in data:
for iter = 1:NrSpikes
    [~,file,~] = fileparts(traceslist{iter});
  	data = load(fullfile(DIR,traceslist{iter}),'traces');
    
    % now that we know the dimensions, setup our matrices
    if iter == 1
        traces_array = zeros(NrSpikes, length(n), size(data.traces,2));
    end
    
    % add data to array
    for ROIiter = 1:length(n);
        traces_array(iter, ROIiter, :) = data.traces(n(ROIiter),:);
    end
    
end

% make difference matrix 'd' & summed dissimilarity vector 'sumd':
for i = 1:NrSpikes
   for j = (i+1):NrSpikes
       % subtract videos from each other ROI by ROI
       for ROIiter = 1:length(n)
           temp = traces_array(i, ROIiter, :) - traces_array(j, ROIiter, :);
           d(i, j) = mean(abs(temp(:))); % MAE: mean absolute error
           traces_sumd(end+1) = d(i, j);
       end
    end
end

% save summed dissimilarity vector for all sleep cells:
save('Dissimilarity_traces','traces_sumd','traces_array','-v7.3')

% plot summed dissimilarity distribution:
figure(2); hold on; title('Dissimilarity histogram traces');
hist(traces_sumd, 200); hold off
saveas(2,'Dissimilarity Histogram traces.eps')

%%

% %% dissimilarity controls
% Non-ROI area in video corrected for area size
% Subtract traces from different cells from each other 
% (should result in larger dissimilarity scores)