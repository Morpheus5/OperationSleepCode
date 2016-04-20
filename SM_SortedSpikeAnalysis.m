close all; 
% Make sure that Spike_Data is loaded (Spontaneous_Data.mat from Byron's
% scripts) for calculation of burst starting & end times

%% Variables:
num_skipped = 1; 

WindowParameter = 4; % in frames 
% so WP * framerate == window in seconds
% Wills data is 22 fr/sec; mine 30 fr/sec

TakeOutLowNrs = 0; % For plotting ROIindex
% Plotted are all connections that fired together more often than this
% value

%% Select which ROIs to plot (ROIindex will still have all ROIs)

num_ROIs = size(Spike_Data,2);

% For Wills singing data:
% n = sort([  2  4 5 7 11 13 14 15 18 19 20 21  24  27 28 30 31  1 32 33]); 
%      1  2 3 4  5  6  7  8  9 10 11 12  13  14 15 16 17 18 19 20

% For SanneSleep6 as comparison to Wills: (14 ROIs / 8&9, 15&16 overlapped)
n = [ 1 2 3 4 5 6 7 8 10 11 12 13 14 15 ];

% For SanneSleep7 all ROIs: (28 ROIs / top left group / 10&11 overlapped)
% n = [ 1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29];

% For SanneSleep7 as comparison to Wills: (21 ROIs / top left group / 10&11 overlapped)
% n = [ 2 3 4 5 6 7 8 9 11 12 13 14 15 18 19 21 22 24 25 26 27];

% For SanneSleep14 as comparison to Wills: (20 ROIs same as Wills 33=29)
% n = sort([  2 4 5 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 30 6 1 3 9]);
% %      1 2 3 4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20

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
% n = 1:num_ROIs;


num_videos = size(Spike_Data,1);
% num_videos = 18

num_frames = size(Spike_Data,3);
num_ROIs = size(n,2);

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

savefilename1 = ['AllVideoSpikingFrames_' num2str(WindowParameter) 'frameswindow'];
save(savefilename1,'AllVideoSpikingFrames')

% AllVideoSpikingFrames shows per video, per ROI, in which frame(s) the cells spiked 
% # values in SpikeData > 0; 


%% Which ROIs spike at the same time? = Find co-activity

SpikeTogetherMatrix = [];
Allspikes = {};

for VideoIter2 = 1:num_videos;

    for FrameIter1 = 1:num_frames;
        [I,J] = find(AllVideoSpikingFrames{VideoIter2} >= FrameIter1 & AllVideoSpikingFrames{VideoIter2} <= FrameIter1+WindowParameter);
        % I = ROI identity that was active in the window
        % Because of the window, it will only pull out ROIs that were
        % active at the same time. 

        % J = which column was it in? to look up the exact frame 
        %       each J corresponds to the I in the same row

        if length(I) > 1;
            num_I = length(I);

            for I_Iter = 1: (num_I) - 1
                frameI = AllVideoSpikingFrames{VideoIter2}(I(I_Iter),J(I_Iter));
                frameJ = AllVideoSpikingFrames{VideoIter2}(I(I_Iter+1),J(I_Iter+1));
                if frameI~=0; 
                    SpikeTogetherMatrix(end+1,1) = I(I_Iter);
                    SpikeTogetherMatrix(end,2) = I(I_Iter+1);
                    SpikeTogetherMatrix(end,3) = (VideoIter2 - 1) * num_frames + frameI;
                    SpikeTogetherMatrix(end,4) = (VideoIter2 - 1) * num_frames + frameJ;
                end
            end
        end
        clear('I','J');
    end
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
    ROIindex(ROIX,ROIY) = ROIindex(ROIX,ROIY)+1;
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

% The first dimension is ROIs (x)
% The second is the ROIs they fire with in the window selected (y)
% The third dimension has the mean lag (z1) and standard deviation (z2)
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

% WARNING: Contains mistake -- The sorted ROIindex matrix is NOT the same 
% as the original one (error might be in sorting or in naming)

AllSeqIndex = [];
ROIseqIndex = [];

% find which ROIs spiked first in each video
for VideoSeqIter = 1:num_videos
    AllVideoSpikingFrames{VideoSeqIter}(AllVideoSpikingFrames{VideoSeqIter}==0)=Inf; 
    % I want the zeros to come last (cells that haven't spiked)
    [SequenceMatrix, SeqIndex] = sortrows((AllVideoSpikingFrames{VideoSeqIter}(n,:)),1);
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
    
    temp = unqASI(valueIndex(1),1); % = the ROI that is most often the first to fire
    
    for i=2:length(valueIndex) % find a next ROI in line if this ROI already fired earlier
        if ismember(temp,SeqIndex)==1
            temp = unqASI(valueIndex(i),1);
        end
        if ismember(temp,SeqIndex)==1
            where = ismember(n,SeqIndex)==0;
            [~,b] = find(where==1);
            temp = n(b(1));
        end
    end
         
    SeqIndex(1,SeqIter1) = temp; % an array with which ROIs fired first  
    % most often, which second, etc., ending with ROIs that did not fire
end

% Sort ROIindex using SeqIndex:
SortedROIindex = [];
for SeqIter3 = SeqIndex
    SortedROIindex(end+1,:) = ROIindex(SeqIter3,:);
end


%% Prepare data for plotting

ROIplotting = SortedROIindex; 
ROIplotting(ROIplotting<TakeOutLowNrs) = 0; % select only more active ROI-pairs
% ROIplotting = ROIplotting - min(ROIplotting(:)); % set range between [0, inf)
% ROIplotting = ROIplotting ./ max(ROIplotting(:)) ; % set range between [0, 1]

LagPlotting = squeeze(ROIlags(:,:,1)) + WindowParameter; % Add WP to make sure all values are positive
% LagPlotting = LagPlotting - min(LagPlotting(:)); % set range between [0, inf)
% LagPlotting = LagPlotting ./ max(LagPlotting(:)) ; % set range between [0, 1]
LagPlotting(isnan(LagPlotting))=0; %replace NaNs by zeros

% Sort LagPlotting using SeqIndex:
SortedLagPlotting = [];
for SeqIter3 = SeqIndex
    SortedLagPlotting(end+1,:) = LagPlotting(SeqIter3,:);
end


LagVar = squeeze(ROIlags(:,:,2)) + WindowParameter; % Add WP to make sure all values are positive
% LagVar = LagVar - min(LagVar(:)); % set range between [0, inf)
% LagVar = LagVar ./ max(LagVar(:)) ; % set range between [0, 1]
LagVar(isnan(LagVar))=0; %replace NaNs by zeros

% Sort LagVar using SeqIndex:
SortedLagVar = [];
for SeqIter3 = SeqIndex
    SortedLagVar(end+1,:) = LagVar(SeqIter3,:);
end


% Remove values that are below the TakeOutLowNrs threshold for sparser
% plotting
ROIXlooping = 1;
ROIYlooping = 1;
AdjustedLagPlotting = [];
AdjustedLagVar = [];

for ROIXlooping = 1:num_ROIs
    for ROIYlooping = 1:num_ROIs
        if ROIplotting(ROIXlooping,ROIYlooping) > 0; 
            AdjustedLagPlotting(ROIXlooping,ROIYlooping) = SortedLagPlotting(ROIXlooping,ROIYlooping); %
            AdjustedLagVar(ROIXlooping,ROIYlooping) = SortedLagVar(ROIXlooping,ROIYlooping); %
        else AdjustedLagPlotting(ROIXlooping,ROIYlooping) = 0;
            AdjustedLagVar(ROIXlooping,ROIYlooping) = 0;
        end
       
    end
end

savefilename2 = ['ROIindexAllData_' num2str(WindowParameter) 'frameswindow'];
save(savefilename2, 'ROIindex','SortedROIindex','SeqIndex','SpikeTogetherMatrix','SpikeTogetherUnique','ROIorder','ROIlags','n','NrOfSpikesIn2ROIs','TakeOutLowNrs','WindowParameter','ROIplotting','AdjustedLagPlotting','AdjustedLagVar') 


%% Plot figures

cmapgray = gray(10*(length(ROIplotting))+10);
cmapgray = flipud(cmapgray(1:10:(end-10),:));

cmapcool = cool(10*(length(AdjustedLagPlotting)));
cmapcool = flipud(cmapcool(1:10:(end),:));

NewNodes = zeros(1,length(n));
nodeLabels = {};
for NodeIter = 1:length(SeqIndex);
    CurrSeq = SeqIndex(1,NodeIter);
    NewNodeName = n(1,CurrSeq);
    nodeLabels{1,NodeIter} = num2str(NewNodeName);
    NewNodes(1,NodeIter) = NewNodeName;
end
 

% Plot & save 1
CoActFigLag = figure('Name','Co-activity index indicated in line width; leading/lagging indicated in grayscale'); 
circularGraph(ROIplotting,'ColorMap',cmapgray,'ColorValue',AdjustedLagPlotting,'Label',nodeLabels);  
colorbar('Position',[0.85 0.15 0.02 0.7],'ColorMap',cmapgray,'Ticks',[0 0.5 1], 'Ticklabels',{'lagging','simultaneous','leading'})

savefilename3 = ['ROI_Co-activity_&_Lag_' num2str(WindowParameter) 'frameswindow'];
saveas(CoActFigLag,savefilename3,'png');
saveas(CoActFigLag,savefilename3,'fig');

% Plot & save 2
CoActFig = figure('Name','Co-activity index indicated in line width and grayscale'); 
circularGraph(ROIplotting,'ColorMap',cmapgray,'ColorValue',ROIplotting,'Label',nodeLabels); 
colorbar('Position',[0.85 0.15 0.02 0.7],'ColorMap',cmapgray,'Ticks',[0 0.5 1], 'Ticklabels',{'low co-activity','medium co-activity','high co-activity'})

savefilename4 = ['ROI_Co-activityIndex_' num2str(WindowParameter) 'frameswindow'];
saveas(CoActFig,savefilename4,'png');
saveas(CoActFig,savefilename4,'fig');

% Plot & save 3
LagFig = figure('Name','Leading/lagging indicated in colormap; Variation in line width'); 
circularGraph(AdjustedLagVar,'ColorMap',cmapcool,'ColorValue',AdjustedLagPlotting,'Label',nodeLabels);
colorbar('Position',[0.85 0.15 0.02 0.7],'ColorMap',cmapcool, 'Ticks', [0 0.5 1], 'Ticklabels',{'lagging','simultaneous','leading'})

savefilename5 = ['ROI_Lags_&_Variation_' num2str(WindowParameter) 'frameswindow'];
saveas(LagFig,savefilename5,'png');
saveas(LagFig,savefilename5,'fig');


% Plot ROI map
 
figure('Name','How often these cells fire together with other cells'); image(ROI.reference_image);
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
savefilename6 = ['ROI_Co-activity_BrainMap_' num2str(WindowParameter) 'frameswindow'];
saveas(ImFig,savefilename6,'png');

%%
autoArrangeFigures(2,3)


