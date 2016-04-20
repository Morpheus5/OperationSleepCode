%% Variables:

num_threshold = 4; % Try 4 to begin with... 
% The smaller this number, the more liberally bursts will be selected
num_skipped = 3; % Try 3 to begin with... 
% The larger this number, the more liberally bursts will be selected
    % if n = 1; it finds the first frame w/o activity
    % if n > 1; it will create an array with the first n frames w/o act.
    % Currently, the script here will select the last frame of the array,
    % so that the burst continues until the n-th time w/o activity. 
num_burstlength = 150; % How long are the longest bursts? I assume none is 
% longer than 2 seconds, 60 frames
num_beginning = 0; % How many frames to include before bursts?

% - Make sure that Av_Data is loaded (Spontaneous_Data.mat from Byron's
% scripts) for calculation of burst starting & end times;
% - Make sure that ROI (roi_data_image.mat) from Will's scripts is loaded 
% for plotting the figure

%% Optionally: Loop through multiple Av_Data files in Processed_Data made by Ca_im_DataStream


% 
% for RecordingIter = 1:size(Processed_Data.Av_Data,2);
%     Av_Data = Processed_Data.Av_Data{1,RecordingIter};
%     
    
    %% End optional part for Processed_Data
    

    %% Looping through videos:

    num_videos = size(Av_Data,1);
    num_frames = size(Av_Data,2);

    for VideoIter = 1: num_videos % loop through videos


    %% PART I: FIND BURSTS
    
    
    [activitylocations] = find(Av_Data(VideoIter,:)>(num_threshold-1)); 
    % find where there is activity of a defined # ROIs simultaneously
    num_bursts = length(activitylocations);
    BurstingFrames = [];
    BurstStartingFrames = [];
    BurstEndingFrames = []; 
    % allocate matrices
  
    % Find which activity frames are individual bursts (eliminate
    % overlapping numbers)
    
    [diff_al] = diff(activitylocations);
    [burstlocations(1,1)] = activitylocations(1,1);
    for i = 2:(length(diff_al))
    if (diff_al(1,i-1)) <= (num_skipped+1); continue; 
    else [burstlocations(end+1)] = activitylocations(1,i);
    end
    end
    num_bursts = length(burstlocations);
    
    
    %% Find beginning of bursts 
  
    for FindBeginningIter = 1:num_bursts % loop through above found bursts
        CurrentBurst1 = burstlocations(1,FindBeginningIter);
        
        Beginning = CurrentBurst1 - num_skipped;
        if Beginning <= 1; Beginning = 1; end
        
        Ending = CurrentBurst1 + num_skipped;
        if Ending >= num_frames; 
            Ending = num_frames; 
            Beginning = num_frames - (2*num_skipped); 
        end
        
        if Beginning == 1; Ending = 2*num_skipped; end
        
        [BurstingFrames] = Av_Data(VideoIter,(Beginning:Ending));
        [~,Y] = find(BurstingFrames,1);
        if size(Y) == [1,0]; Y=0; end;
        
        Beginning_Value = Beginning+Y-num_beginning;
        if Beginning_Value <= 1; Beginning_Value =1; end;
        BurstStartingFrames(end+1) = Beginning_Value;

    end


    % Find ending of bursts 
    
    

    for FindEndingIter = 1:num_bursts 
        if num_bursts <= 1;
        fprintf('\n  \n No bursts survived your criteria in video %0.0f  \n \n',VideoIter);
        end

        CurrentBurst2 = BurstStartingFrames(1,FindEndingIter);
        [~,Y] = find(Av_Data(VideoIter,(CurrentBurst2:end))<1,(num_beginning+num_skipped));
        
        if size(Y) == [1,0]; Y=0; end;
        BurstEndingFrames(end+1) = CurrentBurst2+Y(end);

    end

    % Total data: rows = nr of bursts; column 1 = start, column 2 = end

    AllBurstingFrames = [];
    AllBurstingFrames(1,1)=BurstStartingFrames(1,1);
    AllBurstingFrames(1,2)=BurstEndingFrames(1,1);

    for BurstIter = 2:num_bursts
        
        AllBurstingFrames(BurstIter,2)=BurstEndingFrames(1,BurstIter);
        
        if BurstStartingFrames(1,BurstIter-1) ~= BurstStartingFrames(1,BurstIter)-1;
            AllBurstingFrames(BurstIter,1)=BurstStartingFrames(1,BurstIter);
        else BurstIter = BurstIter+1; end
               
    end

    AllVideoBurstingFrames{VideoIter}=(AllBurstingFrames);
    savefilename = ['AllVideoBurstingFrames_thr' num2str(num_threshold)];%['AllVideoBurstingFrames_rec' num2str(RecordingIter) '_thr' num2str(num_threshold)];
    save(savefilename,'AllVideoBurstingFrames')
      
    clear ('AllBurstingFrames' ,'BurstingFrames' , 'BurstStartingFrames' , 'BurstEndingFrames' , 'Beginning', 'Ending' , 'activitylocations', 'burstlocations', 'num_bursts' , 'diff_al' , 'FindBeginningIter' , 'FindEndingIter' , 'CurrentBurst1' , 'CurrentBurst2' , 'BurstIter' );
    
end % end VideoIter

%end % end RecordingIter


%% PART II: Find ROIs that were involved in bursts

ROIbursting = {};
      
for VideoIter2 = 1:length(AllVideoBurstingFrames);
%   ROItimes = squeeze(Processed_Data.Spike_Data(VideoIter2,:,:));
    ROItimes = squeeze(Spike_Data(VideoIter2,:,:));
        
    for BurstIter2 = 1:size(AllVideoBurstingFrames{1,VideoIter2},1) 
        % Calculates how many bursts there are for each video again
        A = (AllVideoBurstingFrames{1,VideoIter2}(BurstIter2,1));
        B = (AllVideoBurstingFrames{1,VideoIter2}(BurstIter2,2));
        if B >= size(AllVideoBurstingFrames{1,VideoIter2},1); 
            B = size(AllVideoBurstingFrames{1,VideoIter2},1);
        end
        ROIbursting{1,end+1} = ROItimes(:,(A:B)); % saves ROI data for each 
        % burst in the dataset, without indicating from which video it was
    end

% Use the below code lines if you want to save data for each video
% separately (now saves all videos together below the "end" statement):
%     savefilename = ['ROIbursts_Video' num2str(VideoIter2)];
%     save(savefilename,'ROIbursting');
end

save('ROIbursting','ROIbursting')



%% PART III: Analyze which ROIs initiated bursts

BurstParticipation = zeros((size(ROItimes,1)),num_burstlength); 
% Rows: All ROIs 
% Columns: 1=total count of burst participation in frame 1, 2=frame 2, etc


for BurstIter3 = 1:size(ROIbursting,2) % loop through data per burst
    data = ROIbursting{1,BurstIter3};
    num_frames2 = size(data,2);
    
    for ROIiter = 1:size(data,1) % loop through ROIs
        for FrameIter = 1:num_frames2; % loop through each frame in the burst
            if data(ROIiter,FrameIter) == 1
                BurstParticipation(ROIiter,FrameIter) = BurstParticipation(ROIiter,FrameIter)+1;
            elseif data ~= 1 
                FrameIter = FrameIter+1;
            end
        end
    end
end

% Do the same but now calculate the ratio (nr devided by nr of bursts):

BurstParticipationRatio = zeros((size(ROItimes,1)),30); 

data = BurstParticipation;
    
for ROIiter2 = 1:size(data,1) % loop through ROIs
    for FrameIter2 = 1:size(data,2); % loop through each frame in the burst
        if data(ROIiter2,FrameIter2) ~= 0
            BurstParticipationRatio(ROIiter2,FrameIter2) = (data(ROIiter2,FrameIter2))/(size(ROIbursting,2));
        elseif data == 0
            FrameIter2 = (FrameIter2)+1;
        end
    end
end

save('BurstParticipation','BurstParticipation')
save('BurstParticipationRatio','BurstParticipationRatio')


TrueBurstParticipation = sum(ROItimes');
TrueBurstParticipation = TrueBurstParticipation';
save('TrueBurstParticipation','TrueBurstParticipation')

%% Calculate intervals

Intervals = [];

for VideoIter3 = 1:length(AllVideoBurstingFrames);
    for IntervalIter = 1:size(AllVideoBurstingFrames{VideoIter3},1)-1;
        [Start] = AllVideoBurstingFrames{VideoIter3}(IntervalIter,2);
        [End] = AllVideoBurstingFrames{VideoIter3}(IntervalIter+1,1);
        if End-Start <= 1; continue; end;
        Intervals(1,end+1) = End-Start;      
    end
end

MeanIntervals = mean(Intervals);
MedianIntervals = median(Intervals);
StDevIntervals = std(Intervals);
[N,EDGES] = histcounts(Intervals);
% [N,EDGES] = histcounts(X) partitions the values in X into bins, and
    % returns the count in each bin, as well as the bin edges. histcounts
    % determines the bin edges using an automatic binning algorithm that
    % returns uniform bins of a width that is chosen to cover the range of
    % values in X and reveal the shape of the underlying distribution.
    
save('AllVideo_Inter-Burst-Intervals','Intervals')
save('Stats_Inter-Burst-Intervals','MeanIntervals','MedianIntervals','StDevIntervals','N','EDGES')


%% Calculate bursting durations

BurstDurations = [];

for VideoIter4 = 1:length(AllVideoBurstingFrames);
    for IntervalIter = 1:size(AllVideoBurstingFrames{VideoIter4},1)-1;
        [Start] = AllVideoBurstingFrames{VideoIter4}(IntervalIter,1);
        [End] = AllVideoBurstingFrames{VideoIter4}(IntervalIter,2);
        if End-Start <= 1; continue; end;
        BurstDurations(1,end+1) = End-Start;      
    end
end

MeanBurstDurations = mean(BurstDurations);
MedianBurstDurations = median(BurstDurations);
StDevBurstDurations = std(BurstDurations);
[N,EDGES] = histcounts(BurstDurations);
% [N,EDGES] = histcounts(X) partitions the values in X into bins, and
    % returns the count in each bin, as well as the bin edges. histcounts
    % determines the bin edges using an automatic binning algorithm that
    % returns uniform bins of a width that is chosen to cover the range of
    % values in X and reveal the shape of the underlying distribution.
    
save('AllVideo_BurstDurations','BurstDurations')
save('Stats_BurstDurations','MeanBurstDurations','MedianBurstDurations','StDevBurstDurations','N','EDGES')




%% Lastly: Plot figure to visualize which ROIs initiate bursting

% image(ROI.reference_image);
% colormap('gray(2555');
% circle_colors = hot(100);
% hold on;
% 
% NrOfROIs = length(ROI.coordinates);
% for ROIiter3 = 1:NrOfROIs % loop through ROIs
% 
%     %RelativeValues = mean(BurstParticipationRatio,2); % take numbers from BurstParticipationRatios calculated earlier
%     RelativeValues = TrueBurstParticipation;
%     RelativeValues = RelativeValues - min(RelativeValues(:)); % set range between [0, inf)
%     RelativeValues = RelativeValues ./ max(RelativeValues(:)) ; % set range between [0, 1]
%     RelativeValues(RelativeValues==0) = 0.01; % replace 0 by small number
%     ColorValue = RelativeValues(ROIiter3,1);
%     
%     coords = ROI.stats(ROIiter3).Centroid;
%     radius = ROI.stats(ROIiter3).Diameter / 2;
%     scatter(coords(1), coords(2), pi * radius ^ 2, circle_colors(round(ColorValue * 100), :), 'filled');
%         
% % alternative plotting method:
% % top_left = coords - radius; ellipse_size = radius .* [2 2]; 
% % imellipse(gca, [top_left ellipse_size]); % an elipse with draggable edges
% 
% end
