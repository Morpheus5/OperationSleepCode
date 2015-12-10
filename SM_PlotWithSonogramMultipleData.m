%% VARIABLES
% load saved variables into the workspace first:
    % load roi_ave 
    % and the output matrix from plot one or multiple videos 
    % and the stimulus: SONG = audioread('filename');
    
sonogram_comptime = 150; % frame at which song started playing 
comptime = 10; % frames clipped from beginning and compensated in graph
data_to_plot = total_dffvideos; % matrix from which data is grabbed
fs = 48000;
SelectedVideos = [11:29]; 
shiftupValue = 0.8;

startTime = (1*(1/framerate));
endTime = (Frames*(1/framerate));

%% Plot sonogram
[b,a] = ellip(5,.2,80,[500]/(fs/2),'high');
[IMAGE,F,T] = fb_pretty_sonogram(filtfilt(b,a,SONG./abs(max(SONG))),fs,'low',1.5,'zeropad',0);
figure()

    h(1) = subplot(10,1,1:2);
        imagesc((sonogram_comptime+T*framerate-startTime*framerate+1)/framerate,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
        colormap(flipud(bone)*.8)
        freezeColors;
        ylim([0 9000]);
        hold on

    h(2) = subplot(10,1,3:10);
    hold on

%% Plot ROI data  

for videoIter = SelectedVideos 
    plotPerVideo = data_to_plot(:,:,videoIter); 
    
    counter = 1;               
    c = colormap(lines(50)); 
      
        for ROIiter = 1:ROIS
        dataPerROI = plotPerVideo(ROIiter,:); 
        shiftup = (shiftupValue*counter-1); 
        plot(((((comptime+(1:length(dataPerROI)))))/framerate),dataPerROI+shiftup,'Color',c(ROIiter,:),'LineWidth',1.0);
        counter=counter+1;

        end
end 

%% Some more figure stuff
ylabel('');
xlabel('Time (s)');
set(findall(gcf,'type','text'),'FontSize',13,'fontWeight','bold')
linkaxes(h,'x'); 
xlim([0 ((comptime+BOStime)/framerate)]);

