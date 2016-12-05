function [AllVideosBurstingFrames, AllBurstingDurations, ROIbursting, BurstParticipation, BurstParticipationRatio, Intervals, IntervalStats, BurstDurations, DurationStats, ROI_threshold, num_skipped, savefilename_AllVideosBursting, BurstFigure, BurstingImFigure] = SM_BurstStartEndtimes(Av_Data,Spike_Data,ROI,ROI_threshold,num_skipped,num_beginning,min_burstduration)

%% default parameters
if ~exist('ROI_threshold', 'var') || isempty(ROI_threshold)
    ROI_threshold = 4;
end
% The smaller this number, the more liberally bursts will be selected

if ~exist('num_skipped', 'var') || isempty(num_skipped)
    num_skipped = 4;
end
% The larger this number, the more liberally bursts will be selected

if ~exist('num_beginning', 'var') || isempty(num_beginning)
    num_beginning = 0;
end
% How many frames to include before bursts?

if ~exist('min_burstduration', 'var') || isempty(min_burstduration)
    min_burstduration = 2;
end
% As a default, is does not count bursts with length of 1 frame




%% Looping through videos:

num_videos = size(Av_Data,1);
num_frames = size(Av_Data,2);

for VideoIter = 1: num_videos % loop through videos


    %% FIND BURSTS
    
    [activitylocations] = find(Av_Data(VideoIter,:)>=(ROI_threshold)); 
    % find where there is activity of a defined # ROIs simultaneously
    BurstStartingFrames = [];
    BurstEndingFrames = []; 
    % allocate matrices & make sure there's no old data
  
    % Find which activity frames are individual bursts (eliminate
    % overlapping numbers)
    [diff_al] = diff(activitylocations);
    if size(activitylocations) >= [1,1];
        [burstlocations(1,1)] = activitylocations(1,1);
        for i = 2:(length(diff_al))
        if (diff_al(1,i-1)) <= (num_skipped+1); continue; 
        else [burstlocations(end+1)] = activitylocations(1,i);
        end
        end
    else continue;
    end
    num_bursts = length(burstlocations);
    
    
    %% Find beginning of bursts 
  
    for BurstIter = 1:num_bursts % loop through above found bursts
        EarliestActivity = burstlocations(1,BurstIter);
        
        while true
            SearchForActivityStart = EarliestActivity - num_skipped;
            SearchForActivityStop = EarliestActivity - 1;
            
            % stop once we have earliest activity at frame 1
            if SearchForActivityStop < 1
                break;
            end
            
            % do not look before beginning of video
            if SearchForActivityStart < 1
                SearchForActivityStart = 1;
            end
            
            % find earliest during search window
            BurstingFrame = find(Av_Data(VideoIter, SearchForActivityStart:SearchForActivityStop), 1);
            
            % nothing found? stop the loop
            if isempty(BurstingFrame)
                break;
            end
            
            % found something
            EarliestActivity = SearchForActivityStart + BurstingFrame - 1;
        end
        
        BurstStartingFrames(end+1) = EarliestActivity;
    end
    
    % adjust beginnings
    BurstStartingFrames = max(BurstStartingFrames - num_beginning, 1);


    % Find ending of bursts 
  
    for BurstIter = 1:num_bursts % loop through above found bursts
        LatestActivity = burstlocations(1, BurstIter);
        
        while true
            SearchForActivityStart = LatestActivity + 1;
            SearchForActivityStop = LatestActivity + num_skipped;
            
            % stop once we have earliest activity at frame 1
            if SearchForActivityStart > num_frames
                break;
            end
            
            % do not look before beginning of video
            if SearchForActivityStop > num_frames
                SearchForActivityStop = num_frames;
            end
            
            % find earliest during search window
            BurstingFrame = find(Av_Data(VideoIter, SearchForActivityStart:SearchForActivityStop), 1);
            
            % nothing found? stop the loop
            if isempty(BurstingFrame)
                break;
            end
            
            % found something
            LatestActivity = SearchForActivityStart + BurstingFrame - 1;
        end
        
        BurstEndingFrames(end+1) = LatestActivity;
    end

    % Total data: rows = nr of bursts; column 1 = start, column 2 = end

    AllBurstingFrames = [BurstStartingFrames(1) BurstEndingFrames(1)];

    for BurstIter = 2:num_bursts
        if BurstStartingFrames(BurstIter) <= BurstEndingFrames(BurstIter - 1)
            AllBurstingFrames(end, 2) = BurstEndingFrames(BurstIter);
        else
            AllBurstingFrames(end + 1, 1) = BurstStartingFrames(BurstIter);
            AllBurstingFrames(end, 2) = BurstEndingFrames(BurstIter);
        end
    end
    AllBurstingDurations = AllBurstingFrames(:, 2) - AllBurstingFrames(:, 1) + 1;
    AllBurstingFrames = AllBurstingFrames(AllBurstingDurations >= min_burstduration, :);
    AllBurstingDurations = AllBurstingDurations(AllBurstingDurations >= min_burstduration);

    AllVideosBurstingFrames{VideoIter}=(AllBurstingFrames);
    
    savefilename_AllVideosBursting = ['AllVideosBursting_' num2str(ROI_threshold) 'ROIthresh_' num2str(num_skipped) 'FrameWindow' ];
    save(savefilename_AllVideosBursting,'AllVideosBurstingFrames', 'AllBurstingDurations')
      
    clear ( 'activitylocations', 'burstlocations', 'diff_al' );
    
    savefilename = ['BurstFigure_' num2str(ROI_threshold) 'ROIthresh_' num2str(num_skipped) 'FrameWindow' ];
    BurstFigure = figure('Name',['Bursts_Video' num2str(VideoIter)]);
    t = (1:size(Av_Data, 2)) ./ 30; %times
    plot(t, Av_Data); hold on; scatter((BurstStartingFrames(:)./30), (Av_Data(BurstStartingFrames)./30), 'g'); scatter((BurstEndingFrames(:)./30), (Av_Data(BurstEndingFrames)./30), 'r'); hold off;
    saveas(BurstFigure,savefilename,'png');
    saveas(BurstFigure,savefilename,'fig');
     
end % end VideoIter


%% Find ROIs that were involved in bursts

ROIbursting = {};
      
for VideoIter = 1:length(AllVideosBurstingFrames);
    ROItimes = squeeze(Spike_Data(VideoIter,:,:));
        
    for BurstIter = 1:size(AllVideosBurstingFrames{1,VideoIter},1) 
        % Calculates how many bursts there are for each video again
        A = (AllVideosBurstingFrames{1,VideoIter}(BurstIter,1));
        B = (AllVideosBurstingFrames{1,VideoIter}(BurstIter,2));
        ROIbursting{1,end+1} = ROItimes(:,(A:B)); % saves ROI data for each 
        % burst in the dataset, without indicating from which video it was
    end

% Use the below code lines if you want to save data for each video
% separately (now saves all videos together below the "end" statement):
%     savefilename = ['ROIbursts_Video' num2str(VideoIter2)];
%     save(savefilename,'ROIbursting');
end


%% Analyze which ROIs initiated bursts

BurstParticipation = zeros((size(ROItimes,1)), max(AllBurstingDurations)); 
% Rows: All ROIs 
% Columns: 1=total count of burst participation in frame 1, 2=frame 2, etc


for BurstIter = 1:size(ROIbursting,2) % loop through data per burst
    data = ROIbursting{1,BurstIter};
    num_frames2 = size(data,2);
    BurstParticipation(:, 1:num_frames2) = BurstParticipation(:, 1:num_frames2) + data;
end

% Do the same but now calculate the ratio (nr devided by nr of bursts):

BurstParticipationRatio = BurstParticipation;
for BurstIter = 1:size(BurstParticipationRatio, 2)
    BurstParticipationRatio(:, BurstIter) = BurstParticipationRatio(:, BurstIter) ./ sum(AllBurstingDurations >= BurstIter);
end


%% Calculate intervals

Intervals = [];

for VideoIter = 1:length(AllVideosBurstingFrames);
    for IntervalIter = 1:size(AllVideosBurstingFrames{VideoIter},1)-1;
        [Start] = AllVideosBurstingFrames{VideoIter}(IntervalIter,2);
        [End] = AllVideosBurstingFrames{VideoIter}(IntervalIter+1,1);
        if End-Start <= 1; continue; end;
        Intervals(1,end+1) = End-Start;      
    end
end

IntervalStats.MeanIntervals = mean(Intervals);
IntervalStats.MedianIntervals = median(Intervals);
IntervalStats.StDevIntervals = std(Intervals);
[N,EDGES] = histcounts(Intervals);
IntervalStats.N = N;
IntervalStats.EDGES = EDGES;
% [N,EDGES] = histcounts(X) partitions the values in X into bins, and
    % returns the count in each bin, as well as the bin edges. histcounts
    % determines the bin edges using an automatic binning algorithm that
    % returns uniform bins of a width that is chosen to cover the range of
    % values in X and reveal the shape of the underlying distribution.
    

%% Calculate bursting durations

BurstDurations = [];

for VideoIter = 1:length(AllVideosBurstingFrames);
    for IntervalIter = 1:size(AllVideosBurstingFrames{VideoIter},1);
        [Start] = AllVideosBurstingFrames{VideoIter}(IntervalIter,1);
        [End] = AllVideosBurstingFrames{VideoIter}(IntervalIter,2);
        if End-Start <= 1; continue; end;
        BurstDurations(1,end+1) = End-Start;      
    end
end

DurationStats.MeanBurstDurations = mean(BurstDurations);
DurationStats.MedianBurstDurations = median(BurstDurations);
DurationStats.StDevBurstDurations = std(BurstDurations);
[N,EDGES] = histcounts(BurstDurations);
DurationStats.N = N;
DurationStats.EDGES = EDGES;
% [N,EDGES] = histcounts(X) partitions the values in X into bins, and
    % returns the count in each bin, as well as the bin edges. histcounts
    % determines the bin edges using an automatic binning algorithm that
    % returns uniform bins of a width that is chosen to cover the range of
    % values in X and reveal the shape of the underlying distribution.
    
%% Save all values calculated in "AllBurstingInfo":
save('AllBurstingInfo', 'ROIbursting', 'BurstParticipation', 'BurstParticipationRatio', 'Intervals', 'IntervalStats', 'BurstDurations', 'DurationStats')



%% Lastly: Plot figure to visualize which ROIs initiate bursting

BurstingImFigure = figure('Name','ROIs that initiate bursting'); 

image(ROI.reference_image);
colormap('gray(2555');
circle_colors = hot(100);
hold on;

NrOfROIs = length(ROI.coordinates);
for ROIiter = 1:NrOfROIs % loop through ROIs

    RelativeValues = BurstParticipationRatio(:, 1); % take numbers from BurstParticipationRatios calculated earlier
    RelativeValues = RelativeValues - min(RelativeValues(:)); % set range between [0, inf)
    RelativeValues = RelativeValues ./ max(RelativeValues(:)) ; % set range between [0, 1]
    RelativeValues(RelativeValues==0) = 0.01; % replace 0 by small number
    ColorValue = RelativeValues(ROIiter,1);
    
    coords = ROI.stats(ROIiter).Centroid;
    radius = ROI.stats(ROIiter).Diameter / 2;
    scatter(coords(1), coords(2), pi * radius ^ 2, circle_colors(round(ColorValue * 100), :), 'filled');
    text(coords(1) + 2, coords(2) - 7, num2str(ROIiter));
        
end
savefilename = ['ROI_BurstInitiation_BrainMap_' num2str(ROI_threshold) 'ROIthresh_' num2str(num_skipped) 'FrameWindow' ];
saveas(BurstingImFigure,savefilename,'png');
autoArrangeFigures(2,2)

