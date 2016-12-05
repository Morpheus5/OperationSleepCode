
num_skipped = 1667; % equivalent to 5 frames in 30 fr/sec data (10000/30*5)
num_beginning = 0;
min_burstduration = 667; % equivalent to 2 frames in 30 fr/sec data (2/30 *10000)

AllVideosBurstingFrames = {};
%% Looping through videos:

num_videos = size(SM_Spike_Data,1);
num_frames = size(SM_Spike_Data,3);

for VideoIter = 1: num_videos % loop through videos


    %% FIND BURSTS
    
    [activitylocations] = find(SM_Spike_Data(VideoIter,1,:)); 
    % find where there is activity as detected by SM_SpikeDetectorEphys
    BurstStartingFrames = [];
    BurstEndingFrames = []; 
    % allocate matrices & make sure there's no old data
  
    % Find which activity frames are individual bursts (eliminate
    % overlapping numbers)
    [diff_al] = diff(activitylocations);
    if size(activitylocations) >= [1,1];
        [burstlocations(1,1)] = activitylocations(1,1);
        for i = 2:(length(diff_al))
        if (diff_al(i-1,1)) <= (num_skipped+1); 
            continue; 
        else [burstlocations(end+1)] = activitylocations(i,1);
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
            BurstingFrame = find(SM_Spike_Data(1,1, SearchForActivityStart:SearchForActivityStop), 1);
            
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
            BurstingFrame = find(SM_Spike_Data(1,1, SearchForActivityStart:SearchForActivityStop), 1);
            
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
    
    savefilename_AllVideosBursting = ['AllVideosBursting_' num2str(num_skipped) 'FrameWindow' ];
    save(savefilename_AllVideosBursting,'AllVideosBurstingFrames', 'AllBurstingDurations')
      
    clear ( 'activitylocations', 'burstlocations', 'diff_al' );
    
    savefilename = ['BurstFigure_' num2str(num_skipped) 'FramesSkipped_' num2str(min_burstduration) 'MinDuration' ];
    BurstFigure = figure('Name','Bursts Jeff''s Ephys sleep data' );
    T = 1: 1: size(SM_Spike_Data, 3); % duration in # frames
    plot(T, squeeze(SM_Spike_Data(1,1,:))); 
    hold on;
    plot(T,alldata{1}*10);
    scatter((AllBurstingFrames(:,1)), squeeze(SM_Spike_Data(1,1,AllBurstingFrames(:,1))), 'g'); 
    scatter((AllBurstingFrames(:,2)), squeeze(SM_Spike_Data(1,1,(AllBurstingFrames(:,2)))), 'r'); hold off;
    ylim([-1.5 1.1]); hold off;

    saveas(BurstFigure,savefilename,'png');
    saveas(BurstFigure,savefilename,'fig');
     
end % end VideoIter


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
save('AllBurstingInfo', 'Intervals', 'IntervalStats', 'BurstDurations', 'DurationStats')

