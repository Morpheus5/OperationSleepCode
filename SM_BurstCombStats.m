%% Make BurstingStats folder:
mat_dir = 'CombinedBurstingStats';
if exist(mat_dir,'dir'); 
    rmdir(mat_dir,'s'); 
end
mkdir(mat_dir);
DIR=pwd; 

%% Read in all datafiles:
mov_listing = dir(fullfile(DIR,'Bursting*'));
mov_listing = {mov_listing(:).name}; 

BurstsTotal = [];
IntervalsTotal = [];
  
[nblanks formatstring]=progressbar(length(mov_listing));
fprintf(1,['Progress:  ' blanks(nblanks) ]);

%% Loop through datafiles and concatenate burst durations and intervals:
for videoIter = 1:length(mov_listing);
    
    fprintf(1,[ num2str(videoIter), ' / ', num2str(length(mov_listing)),'\n' ] );
    
    datafile = [ mov_listing{videoIter}, '/', 'AllBurstingInfo.mat' ];
  	load(fullfile(DIR,datafile), 'BurstDurations', 'Intervals' );

    a = length(BurstsTotal);
    b = length(BurstDurations);
    c = length(IntervalsTotal);
    d = length(Intervals);
    
    BurstsTotal(1,a+1:a+b) = BurstDurations(:);
    IntervalsTotal(1,c+1:c+d) = Intervals(:);
end

DurationStats.MeanBurstDurations = mean(BurstsTotal);
DurationStats.MedianBurstDurations = median(BurstsTotal);
DurationStats.StDevBurstDurations = std(BurstsTotal);
[N,EDGES] = histcounts(BurstsTotal);
DurationStats.N = N;
DurationStats.EDGES = EDGES;

IntervalStats.MeanIntervals = mean(IntervalsTotal);
IntervalStats.MedianIntervals = median(IntervalsTotal);
IntervalStats.StDevIntervals = std(IntervalsTotal);
[N,EDGES] = histcounts(IntervalsTotal);
IntervalStats.N = N;
IntervalStats.EDGES = EDGES;   

save('CombinedBurstingStats/Stats','IntervalStats','DurationStats')
save('CombinedBurstingStats/Durations','BurstsTotal')
save('CombinedBurstingStats/Intervals','IntervalsTotal')

edges = [0:0.1:70];
BurstFigure = figure('Name','Burst durations histogram'); 
hold on; 
hist((BurstsTotal(:)./30),edges);
ylim([0 5]); xlim([0 10]);
h = findobj(BurstFigure,'Type','patch');
h.FaceColor = [ 0.6 0.8 0.6 ];
h.EdgeColor = [ 0 0.2 0 ];

edges = [0:0.5:70];
IntervalFigure = figure('Name','Inter-burst-interval durations histogram'); 
hold on; 
hist((IntervalsTotal(:)./30),edges);
ylim([0 5]); xlim([0 40]);
g = findobj(IntervalFigure,'Type','patch');
g.FaceColor = [ 0.4 0.8 0.8 ];
g.EdgeColor = [ 0 0.2 0.2 ];

saveas(BurstFigure,'CombinedBurstingStats/BurstHistogram','png');
saveas(IntervalFigure,'CombinedBurstingStats/IntervalHistogram','png');

