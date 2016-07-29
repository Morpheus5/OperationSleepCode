function [AllVideoSpikingFrames, SpikeTogetherUnique, SortedROIindex, ROIlags, CoActFigLag, CoActFig, LagFig, ImFig] = SM_SortedSpikeAnalysis(Spike_Data,ROI,WP)
% Make sure that Spike_Data is loaded 

WindowParameter = WP; % in frames - try 4 to begin with
% so WP / framerate == window in seconds
% Wills data is 22 fr/sec; mine 30 fr/sec
TakeOutLowNrs = 0; % For plotting ROIindex
% Plotted are all connections that fired together more often than this
% value

%% Select which ROIs to plot (ROIindex will still have all ROIs)

% For Wills singing data:
% n = sort([  2  4 5 7 11 13 14 15 18 19 20 21  24  27 28 30 31  1 32 33]); 

% For SanneSleep6: (14 ROIs / 8&9, 15&16 overlapped)
% n = [ 1 2 3 4 5 6 7 8 10 11 12 13 14 15 ];

% For SanneSleep7 all ROIs: (28 ROIs / top left group / 10&11 overlapped)
% n = [ 1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

% For SanneSleep7 as comparison to Wills: (21 ROIs / top left group / 10&11 overlapped)
% n = [ 2 3 4 5 6 7 8 9 11 12 13 14 15 18 19 21 22 24 25 26 27];

% For SanneSleep14 as comparison to Wills: (20 ROIs same as Wills 33=29)
% n = sort([  2 4 5 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 30 6 1 3 9]);

% For SanneSleep16 (26 ROIs - 21 & 22 overlapped)
% n = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23 24 25 26 27];

% For SanneSleep28 (21 ROIs - 4&5, 8&9, 10&11 overlapped)
% n = [ 1 2 3 4 6 7 8 11 12 13 14 15 16 17 18 19 20 21];

% For SanneSleep31 (1&2, 4&5, 21&22, 34&35 overlapped)
% n = [ 1 3 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27 28 29 30 31 32 33 34 36 ];
% 20 ROIs picked somewhat randomly (skipped about every other ROI):
% n = [ 1 3 5 8 9 10 12 14 16 18 20 23 25 28 29 30 32 33 34 36 ];
% 20 ROIs picked on location ("middle" cluster):
% n = [ 6 7 8 9 16 17 18 19 20 21 25 26 27 28 29 30 31 32 34 36 ];
% n = [ 6 7 8 13 14 15 16 17 18 19 20 21 23 24 25 26 28 29 30 31 ];

% For plotting all ROIs:
num_ROIs = size(Spike_Data,2);
n = 1:num_ROIs;

%% Calculating sizes
num_ROIs = size(n,2);
num_videos = size(Spike_Data,1);
num_frames = size(Spike_Data,3);

%% Find all frames in which spiking occurred

AllSpikingFrames = zeros(num_ROIs,num_frames);
        
for VideoIter1 = 1: num_videos 
    
    for ROIiter1 = n  
        [activitylocations] = find(Spike_Data(VideoIter1,ROIiter1,:)>0); 
        num_spikes = length(activitylocations);
                
            for SpikeIter = 1:num_spikes
                AllSpikingFrames(ROIiter1,SpikeIter)=activitylocations(SpikeIter,1);
            end
            
        clear('activitylocations');
        
    end     
    
    AllVideoSpikingFrames{VideoIter1} = AllSpikingFrames;
    
end

% savefilename1 = ['SpikingFrames_' num2str(WindowParameter) 'frameswindow'];
% save(savefilename1,'AllVideoSpikingFrames')

%% Which ROIs spike at the same time? = Find co-activity

SpikeTogetherMatrix = [];
Allspikes = {};

for video = 1:num_videos
    % rows = ROIs; columns = frames
    spikes = squeeze(Spike_Data(video, :, :));
    max_frames = size(spikes, 2);
    
    % find all spikes
    [rois, frames] = find(spikes);
    
    % for everything that spikes
    for i = 1:length(rois)
        % get the subsequent window
        % (everything that spikes in the currrent frame plus the next WindowParameter - 1 frames)
        subsequent_spikes = spikes(:, frames(i):min(frames(i) + WindowParameter - 1, max_frames));
        
        % prevent double counting (this prevents counting spiking with self
        % in the same frame, and prevents counting cospiking rois in the
        % same frame multiple times)
        subsequent_spikes(1:rois(i), 1) = 0;
        
        % roi and frame; everything else in subsequent spikes is a unique
        % "spike together" with this roi and frame
        roi = rois(i);
        frame = frames(i);
        
        [co_rois, co_frames] = find(subsequent_spikes);
        co_frames = co_frames + frame - 1;
        
        if ~isempty(co_rois)
            new_rows = [ones(length(co_rois), 1) * roi, co_rois, ones(length(co_rois), 1) * frame, co_frames];
            SpikeTogetherMatrix = [SpikeTogetherMatrix; new_rows];
        end
    end
end

if isempty(SpikeTogetherMatrix)
    warning('No spikes detected');
end

I = SpikeTogetherMatrix(:,1);
J = SpikeTogetherMatrix(:,2);
frameI = SpikeTogetherMatrix(:,3);
frameJ = SpikeTogetherMatrix(:,4);
SpikeTogetherTable = table(I,J,frameI,frameJ);
SpikeTogetherUnique = unique(SpikeTogetherTable);

%% Now, make matrices that add up these spikes for each ROI-pair
%  and calculate the lag between the ROI-pairs:

ROIindex = zeros(num_ROIs,num_ROIs);
ROIorder = cell(num_ROIs,num_ROIs);
ROIlags = zeros(num_ROIs,num_ROIs,2);

% To have same nr of spikes as will's singing data & somewhat equally 
% distributed over ROIs, run this instead of the first line below:
% 
% SpikeInterval = floor(size(SpikeTogetherUnique,1)/5688);
% for Xiter = 1 : SpikeInterval : SpikeInterval*5688 

for Xiter = 1 : size(SpikeTogetherUnique,1)
    ROIXcoordinate = SpikeTogetherUnique{Xiter,1};
    [~,ROIX,~] = find(n==ROIXcoordinate); 
    ROIYcoordinate = SpikeTogetherUnique{Xiter,2};
    [~,ROIY,~] = find(n==ROIYcoordinate);
    
%     % OPTION 1: asymmetric; ignore ROIs spiking in same frame
%     if SpikeTogetherUnique{Xiter, 3} == SpikeTogetherUnique{Xiter, 4}
%         continue;
%     end
%     ROIindex(ROIX,ROIY) = ROIindex(ROIX,ROIY)+1;
    
%     % OPTION 2: assymetric; but count spiking together in same frame on
%     % both side
%     if SpikeTogetherUnique{Xiter, 3} == SpikeTogetherUnique{Xiter, 4}
%         ROIindex(ROIY,ROIX) = ROIindex(ROIY,ROIX)+1;
%     end
%     ROIindex(ROIX,ROIY) = ROIindex(ROIX,ROIY)+1;
    
    % OPTION 3: symmetric
    ROIindex(ROIX,ROIY) = ROIindex(ROIX,ROIY)+1;
    ROIindex(ROIY,ROIX) = ROIindex(ROIY,ROIX)+1;
    
    ROIorder{ROIX,ROIY}{end+1,1} = (SpikeTogetherUnique{Xiter,4}) - (SpikeTogetherUnique{Xiter,3});

end

% loop through ROIorder cells (X & Y); 
% then loop through each lag value in those cells 
% (number depends on # co-activities = LagIter)

for XROIorder = 1:num_ROIs
    for YROIorder = 1:num_ROIs

        Lag = {};
        LagMatrix = [];

        for LagIter = 1:length(ROIorder{XROIorder,YROIorder})
            Lag(end+1,1) = (ROIorder{XROIorder,YROIorder}(LagIter,1));
        end

        LagMatrix = cell2mat(Lag);
        ROIlags(XROIorder,YROIorder,1) = mean(LagMatrix);
        ROIlags(XROIorder,YROIorder,2) = std(LagMatrix);

    end
end

NrOfSpikesIn2ROIs = sum(sum(ROIindex(:)));

% ROIlags info: 

% The first dimension (x) is ROIs 
% The second (y) is the ROIs they fire with in the window selected 
% The third dimension (z) has the mean lag (z1) and standard deviation (z2)
    % If there is a positive mean lag with a low std, then cell(x) always 
    % and robustly leads cell(y);
% If there is a mean lag of around zero with a low std, then cell(x) always 
% and robustly fires at the same time as cell(y);
    % If there is a mean lag of around zero with a high std, then cell(x)  
    % might lead sometimes, while it follows cell(y) at other times;
% If there is a negative lag with a low std, then cell(x) always and
% robustly follows cell(y)
    % If the value is NaN (not a number), than these ROIs did not fire 
    % together



%% Stereotyped sequence or random ordering?

% Determine which ROIs spiked first -- the order will be used for sorting 
% the ROIindex as well later (below in the loops)

AllSeqIndex = [];
ROIseqIndex = [];

% find which ROIs spiked first in each video
for VideoSeqIter = 1:num_videos
    AllVideoSpikingFrames{VideoSeqIter}(AllVideoSpikingFrames{VideoSeqIter}==0)=Inf; 
    % I want the zeros to come last (cells that haven't spiked)
    [SequenceMatrix, SeqIndex] = sortrows((AllVideoSpikingFrames{VideoSeqIter}(:,:)),1);
    AllSeqIndex(VideoSeqIter,:) = SeqIndex;
    AllVideoSpikingFrames{VideoSeqIter}(AllVideoSpikingFrames{VideoSeqIter}==Inf)=0; 
end

% find the most common ones if there is variation
SeqIndex = [];
for SeqIter1 = 1:size(AllSeqIndex,2)
    unqASI = unique(AllSeqIndex(:,SeqIter1)); 
    % which ROIs ever fired first, or second, etc?
    [valueCount,~] = histc(AllSeqIndex(:,SeqIter1),unqASI); 
    % and how often did each of them fire?
    [valueSort,valueIndex] = sort(valueCount,'descend'); 
    % sort to find the maximum

    CurrComparison = unqASI(valueIndex(1),1); % = the ROI that is most often the first to fire

    for i=2:length(valueIndex) % find a next ROI in line if this ROI already fired earlier
        if ismember(CurrComparison,SeqIndex)==1
            CurrComparison = unqASI(valueIndex(i),1);
        end
        if ismember(CurrComparison,SeqIndex)==1
            where = ismember(n,SeqIndex)==0;
            [~,b] = find(where==1);
            CurrComparison = n(b(1));
        end
    end

    SeqIndex(1,SeqIter1) = CurrComparison; % an array with which ROIs fired first  
    % most often, which second, etc., ending with ROIs that did not fire
end

%% Sorting based on brain location:
% Find smallest X & smallest Y pair first = top left ROI
% Sort rest of ROIs based on distance from first coordinates

num_n = length(n);

% Make ROI location matrix 'lookuptable':
CoordinateLookup = zeros(num_n,3); % allocate matrix; 

for CoordIter = 1:num_n;
    CoordinateLookup(CoordIter,1) = n(CoordIter); % columns 1-3: ROInr, X, Y
    CoordinateLookup(CoordIter,2:3) = ROI.stats(n(CoordIter)).Centroid;
end

% Make matrix in which ROIs are sorted on distance from the top left ROI:
CoordinateSorting = zeros(num_n,4); % columns: ROInr, X, Y, distance to first ROI
d = zeros(num_n,2); % matrix for distances
LocIndex = zeros(1,num_n);

% Start with top left ROI:
[X,I] = min(CoordinateLookup(:,2)); % = top left ROI
CoordinateSorting(1,1) = n(I(1));
CoordinateSorting(1,2) = CoordinateLookup(I(1),2);
CoordinateSorting(1,3) = CoordinateLookup(I(1),3);
CoordinateSorting(1,4) = 0;
LocIndex(1,1) = I(1);


% Find nearest ROI with precedingly ordered ROI:
for CoordIter2 = 2:num_n;
    Coord1 = [CoordinateSorting(CoordIter2-1,2),CoordinateSorting(CoordIter2-1,3)]; % bijv ROI 16

    % then find distances between all ROIs and previously ordered ROI in d:
    for CoordIter3 = 1:num_n;
        Coord2 = [CoordinateLookup(CoordIter3,2),CoordinateLookup(CoordIter3,3)];
        CurrComparison = [Coord1;Coord2];
        d(CoordIter3,1) = n(CoordIter3);
        d(CoordIter3,2) = pdist(CurrComparison,'euclidean');
    end
    clear('X','I');
    [X,I] = sort(d(:,2)); % nearest ROI - X is distance, n(I) is ROI 

    % find first ROI that is not yet in the Sorting matrix:
    % Find first n(I) isnotmember CoordinateSorting(:,1)
    position = ismember(n(I),CoordinateSorting(:,1));
    position = find(position==0,1);

    % Fill in next ROI in Sorting matrix:
    CoordinateSorting(CoordIter2,1) = n(I(position)); % ROI
    CoordinateSorting(CoordIter2,2) = CoordinateLookup(I(position),2);
    CoordinateSorting(CoordIter2,3) = CoordinateLookup(I(position),3);
    CoordinateSorting(CoordIter2,4) = X(position); % distance

    LocIndex(1,CoordIter2) = I(position);

end

% LocIndex = CoordinateSorting(:,1)';

%% Sort ROIindex using SeqIndex:
%SeqIndex = n; % undo all sorting
SeqIndex = symrcm(ROIindex);

SortedROIindex = ROIindex(SeqIndex, SeqIndex);
LocationROIindex = ROIindex(LocIndex, LocIndex);


%% Sort ROIplotting using SeqIndex:
% ROIplotting = SortedROIindex; 
ROIplotting = LocationROIindex;


%% Prepare data for plotting
ROIplotting(ROIplotting<TakeOutLowNrs) = 0; % select only more active ROI-pairs
% ROIplotting = ROIplotting - min(ROIplotting(:)); % set range between [0, inf)
% ROIplotting = ROIplotting ./ max(ROIplotting(:)) ; % set range between [0, 1]

LagPlotting = squeeze(ROIlags(:,:,1)) + WindowParameter; % Add WP to make sure all values are positive
% LagPlotting = LagPlotting - min(LagPlotting(:)); % set range between [0, inf)
% LagPlotting = LagPlotting ./ max(LagPlotting(:)) ; % set range between [0, 1]
LagPlotting(isnan(LagPlotting))=0; %replace NaNs by zeros


%% Sort LagPlotting using SeqIndex:
LagPlotting = LagPlotting(SeqIndex, SeqIndex);
%LagPlotting = LagPlotting(LocIndex, LocIndex); 

LagVar = squeeze(ROIlags(:,:,2)) + WindowParameter; % Add WP to make sure all values are positive
% LagVar = LagVar - min(LagVar(:)); % set range between [0, inf)
% LagVar = LagVar ./ max(LagVar(:)) ; % set range between [0, 1]
LagVar(isnan(LagVar))=0; %replace NaNs by zeros


%% Sort LagVar using SeqIndex:
% LagVar = LagVar(SeqIndex, SeqIndex);
LagVar = LagVar(LocIndex, LocIndex);

% savefilename2 = ['Video6_ROIindexAllData_' num2str(WindowParameter) 'frameswindow'];
% save(savefilename2, 'ROIindex','SortedROIindex','SeqIndex','SpikeTogetherMatrix','SpikeTogetherUnique','ROIorder','ROIlags','n','NrOfSpikesIn2ROIs','TakeOutLowNrs','WindowParameter','ROIplotting','AdjustedLagPlotting','AdjustedLagVar') 

%% Plot figures

cmapgray = gray(10*(length(ROIplotting))+10);
cmapgray = flipud(cmapgray(1:10:(end-10),:));

cmapcool = cool(10*(length(LagPlotting)));
cmapcool = flipud(cmapcool(1:10:(end),:));

NewNodes = zeros(1,length(n));
nodeLabels = {};
for NodeIter = 1:length(LocIndex);
    CurrSeq = LocIndex(1,NodeIter);
    NewNodeName = n(1,CurrSeq);
    nodeLabels{1,NodeIter} = num2str(NewNodeName);
    NewNodes(1,NodeIter) = NewNodeName;
end


% Plot & save 1
CoActFigLag = figure('Name','Co-activity index indicated in line width; leading/lagging indicated in grayscale'); 
circularGraph(ROIplotting,'ColorMap',cmapgray,'ColorValue',LagPlotting,'ColorRange',[0 2*WindowParameter],'Label',nodeLabels);  
colorbar('Position',[0.85 0.15 0.02 0.7],'ColorMap',cmapgray,'Ticks',[0 0.5 1], 'Ticklabels',{'lagging','simultaneous','leading'})

% savefilename3 = ['ROI_Co-activity_&_Lag_' num2str(WindowParameter) 'frameswindow'];
% saveas(CoActFigLag,savefilename3,'png');
% saveas(CoActFigLag,savefilename3,'fig');

% Plot & save 2
CoActFig = figure('Name','Co-activity index indicated in line width and grayscale'); 
circularGraph(ROIplotting,'ColorMap',cmapgray,'ColorValue',ROIplotting,'Label',nodeLabels); 
colorbar('Position',[0.85 0.15 0.02 0.7],'ColorMap',cmapgray,'Ticks',[0 0.5 1], 'Ticklabels',{'low co-activity','medium co-activity','high co-activity'})

% savefilename4 = ['ROI_Co-activityIndex_' num2str(WindowParameter) 'frameswindow'];
% saveas(CoActFig,savefilename4,'png');
% saveas(CoActFig,savefilename4,'fig');

% Plot & save 3
LagFig = figure('Name','Leading/lagging indicated in colormap; Variation in line width'); 
circularGraph(LagVar,'ColorMap',cmapcool,'ColorValue',LagPlotting,'ColorRange',[0 2*WindowParameter],'Label',nodeLabels);
colorbar('Position',[0.85 0.15 0.02 0.7],'ColorMap',cmapcool, 'Ticks', [0 0.5 1], 'Ticklabels',{'lagging','simultaneous','leading'})
% 
% savefilename5 = ['ROI_Lags_&_Variation_' num2str(WindowParameter) 'frameswindow'];
% saveas(LagFig,savefilename5,'png');
% saveas(LagFig,savefilename5,'fig');


% Plot ROI map

figure('Name','How often these cells fire together with other cells'); 
image(ROI.reference_image);
colormap('gray(2555');
circle_colors = parula(100);
hold on;

NrOfROIs = length(ROI.coordinates);
for ROIfigIter = 1:length(n)

    Subselection = n(1,ROIfigIter); % select which ROIs will be plotted

    Raw = ROIindex;
    TakeOutLowNrs2 = mean(Raw(:));
    Raw(Raw<TakeOutLowNrs2) = 0;
    NewROIindex = mean(Raw);
    RelativeValues = (NewROIindex');
    RelativeValues = RelativeValues - min(RelativeValues(:)); % set range between [0, inf)
    RelativeValues = RelativeValues ./ max(RelativeValues(:)) ; % set range between [0, 1]
    RelativeValues(RelativeValues<=0.01) = 0.01; % replace 0 by small number
    ColorValue = RelativeValues(ROIfigIter,1);

    coords = ROI.stats(Subselection).Centroid;
    radius = ROI.stats(Subselection).Diameter / 2;
    ImFig = scatter(coords(1), coords(2), pi * radius ^ 2, circle_colors(round(ColorValue * 100), :), 'filled');
    text(coords(1) + 2, coords(2) - 7, num2str(Subselection));
end
colorbar('ColorMap',parula)

% savefilename6 = ['ROI_Co-activity_BrainMap_' num2str(WindowParameter) 'frameswindow'];
% saveas(ImFig,savefilename6,'png');

%%
autoArrangeFigures(5,5)
