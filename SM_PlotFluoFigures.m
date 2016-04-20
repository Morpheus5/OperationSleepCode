function SM_PlotFluoFigures(roi_ave)

%% Prepare reading in the data

[Videos]=length(roi_ave.raw);
[Frames]=length(roi_ave.raw{1}(1,:));
[ROIS]=length(roi_ave.raw{1}(:,1));

%% VARIABLES make sure to adjust if necessary

VideosSelected = 1:Videos; %[7 10 13 15 17 ]; 
framerate = 30; % nr of frames per second
compensationtime = 30; % frames to delete at beginning & end (similar to Byron's plots)
BOSselection = compensationtime:(Frames-compensationtime);
PlottedNrOfFrames = length(BOSselection)
savefilename = ['fluoplot_' datestr(clock)];

%% Calculate values and preallocate matrix
BOStime = length(BOSselection);

data = zeros(ROIS,Frames,Videos);
total_raw = zeros(ROIS,BOStime,Videos); 
total_dff = zeros(ROIS,BOStime,Videos); 
total_videos = zeros(ROIS,BOStime,Videos);
total_dffvideos = zeros(ROIS,BOStime,Videos);

%% Read in data
for videoIter = VideosSelected
data(:,:,videoIter)=(roi_ave.raw{videoIter});
end

%% loop through data to calculate DFF

c = colormap(lines(ROIS*Videos+10));

for videoIter2=VideosSelected
    dataPerVideo = data(:,BOSselection,videoIter2);
          
    for ROIiter=1:ROIS
    
    dataPerROI = dataPerVideo(ROIiter,:);
    baseline=prctile(dataPerROI,0); % 0 is  the min
    dff=(dataPerROI-baseline)./baseline;
    total_raw(ROIiter,:,videoIter2)=dataPerROI;
    total_dff(ROIiter,:,videoIter2)=dff;
    
    end
    
    total_videos=total_videos+total_raw;
    total_dffvideos=total_dffvideos+total_dff;
end

%% Plotting 

for videoIter3 = VideosSelected
    plotPerVideo = total_dffvideos(:,:,videoIter3);
    figure(); hold on
    counter = 1;               
     
    for ROIiter2 = 1:ROIS    
        dff = plotPerVideo(ROIiter2,:);
        shiftup = (0.05*counter-1); 
        plot(((0:length(dff))/framerate),dff+shiftup,'Color',c(ROIiter2,:),'LineWidth',1.0);
        xlim([0 (BOStime/framerate)]);
        counter = counter+1;
    end
    hold off

%% More figure stuff & save data and variables
ylabel('');
xlabel('Time (s)');
set(findall(gcf,'type','text'),'FontSize',13,'fontWeight','bold')
  
save(savefilename);
  
end
end