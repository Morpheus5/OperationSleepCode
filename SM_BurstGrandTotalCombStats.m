%% Make BurstingStats folder:
mat_dir = 'GrandTotalStats';
if exist(mat_dir,'dir'); 
    rmdir(mat_dir,'s'); 
end
mkdir(mat_dir);
DIR=pwd; 

%% Read in all datafiles:
mov_listing = dir(fullfile(DIR,'CombinedBurstingStats#*'));
mov_listing = {mov_listing(:).name}; 

TotalDurations = [];
TotalIntervals = [];
  
[nblanks formatstring]=progressbar(length(mov_listing));
fprintf(1,['Progress:  ' blanks(nblanks) ]);

%% Loop through datafiles and concatenate burst durations and intervals:
for videoIter = 1:length(mov_listing);
    
    fprintf(1,[ num2str(videoIter), ' / ', num2str(length(mov_listing)),'\n' ] );
    datafile1 = [ mov_listing{videoIter}, '/'  'Durations'];
  	load(fullfile(DIR,datafile1) );
    datafile2 = [ mov_listing{videoIter}, '/'  'Intervals'];
    load(fullfile(DIR,datafile2) );

    a = length(TotalDurations);
    b = length(BurstsTotal);
    c = length(TotalIntervals);
    d = length(IntervalsTotal);
    TotalDurations(1,a+1:a+b) = BurstsTotal(:);
    TotalIntervals(1,c+1:c+d) = IntervalsTotal(:);
    
end

DurationStats.MeanBurstDurations = mean(TotalDurations);
DurationStats.MedianBurstDurations = median(TotalDurations);
DurationStats.StDevBurstDurations = std(TotalDurations);
[N,EDGES] = histcounts(TotalDurations);
DurationStats.N = N;
DurationStats.EDGES = EDGES;

IntervalStats.MeanIntervals = mean(TotalIntervals);
IntervalStats.MedianIntervals = median(TotalIntervals);
IntervalStats.StDevIntervals = std(TotalIntervals);
[N,EDGES] = histcounts(TotalIntervals);
IntervalStats.N = N;
IntervalStats.EDGES = EDGES;   

save('GrandTotalStats/Stats','IntervalStats','DurationStats')
save('GrandTotalStats/Durations','TotalDurations')
save('GrandTotalStats/Intervals','TotalIntervals')

edges = [0:0.1:70];
BurstFigure = figure('Name','Burst durations histogram'); hold on; hist((IsoDurations(:)./30),edges);
ylim([0 70]); xlim([0 10]);
h = findobj(BurstFigure,'Type','patch');
h.FaceColor = [ 0.6 0.8 0.6 ];
h.EdgeColor = [ 0 0.2 0 ];

edges = [0:0.5:70];
IntervalFigure = figure('Name','Inter-burst-interval durations histogram'); hold on; hist((IsoIntervals(:)./30),edges);
ylim([0 70]); xlim([0 70]);
g = findobj(IntervalFigure,'Type','patch');
g.FaceColor = [ 0.4 0.8 0.8 ];
g.EdgeColor = [ 0 0.2 0.2 ];

saveas(BurstFigure,'GrandTotalStats/BurstHistogram','png');
saveas(IntervalFigure,'GrandTotalStats/IntervalHistogram','png');


