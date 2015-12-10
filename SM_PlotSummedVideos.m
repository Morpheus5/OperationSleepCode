function SM_PlotSummedVideos(roi_ave)

%% prepare reading in data

[Videos]=length(roi_ave.raw);
[Frames]=length(roi_ave.raw{1}(1,:));
[ROIS]=length(roi_ave.raw{1}(:,1));

%% VARIABLES make sure to adjust if necessary
VideosSelected = 1:Videos;  
framerate = 30;
comptime=10; % frames to delete at beginning & compensate in graph
BOSselection = comptime:Frames;
savefilename = ['summed_matrix_' datestr(clock)];

%% Preallocate matrices:

startTime = (1*(1/framerate));
endTime = (Frames*(1/framerate));
BOStime = length(BOSselection);

data=zeros(ROIS,Frames,Videos);
summed_raw=zeros(ROIS,BOStime); 
summed_dff=zeros(ROIS,BOStime); 
summed_videos=zeros(ROIS,BOStime);
summed_dffvideos=zeros(ROIS,BOStime);

%% read in data

for videoIter=VideosSelected
data(:,:,videoIter)=(roi_ave.raw{videoIter});
end

%% loop through data to calculate DFF and make data matrices

for videoIter=VideosSelected 
    
dataPerVideo = data(:,BOSselection,videoIter);

    for ROIiter=1:ROIS
    
        dataPerROI = dataPerVideo(ROIiter,:);
        baseline=prctile(dataPerROI,0); % 0 is  the min
        dff=(dataPerROI-baseline)./baseline;
        summed_raw(ROIiter,:)=dataPerROI;
        summed_dff(ROIiter,:)=dff;
    end
        summed_videos=summed_videos+summed_raw;
        summed_dffvideos=summed_dffvideos+summed_dff;
end


    
%% plot data per ROI    
figure(); hold on
c = colormap(lines(ROIS+10)); 
counter = 1; 
         
    for ROIiter2=1:ROIS    
        dff=summed_dffvideos(ROIiter2,:);
        shiftup = (0.05*counter-1); 
        plot(((comptime+(1:length(dff)))/framerate),dff+shiftup,'Color',c(ROIiter2,:),'LineWidth',1.0);
        xlim([0 ((comptime+BOStime)/framerate)]);
        counter=counter+1;
    end
    
%% More figure stuff
ylabel('');
xlabel('Time (s)');
set(findall(gcf,'type','text'),'FontSize',13,'fontWeight','bold')

%% Save variables and matrices
save(savefilename);
  
end



%% ----[ Setup Sonogram ]------%

% fs=24.4141e3;
% [b,a]=ellip(5,.2,80,[500]/(fs/2),'high');
% [IMAGE,F,T]=fb_pretty_sonogram(filtfilt(b,a,mic_data./abs(max(mic_data))),fs,'low',1.5,'zeropad',0);

%% ----[ Plotting ]------%
% h(1) = subplot(10,1,1:2);
%       imagesc((T*framerate - startTime*framerate + 1)/framerate,F,log(abs(IMAGE)+1e+2));set(gca,'YDir','normal');
%       colormap(flipud(bone)*.8)
%       freezeColors;
%       ylim([0 9000]);
%       hold on

% here, you can make a subplot on top with a sonogram, for example
% h(2) = subplot(10,1,3:10);
% title('Calcium transients LNY19RB 9/2/15 BOS playback #4 (normalized ROIs)');


%% ---------[ Optional: picking peaks ] --------------%

% This code is for selecting peaks in fluorescense change and constructing
% a raster plot based on those
%      
%   [pks,locs] = findpeaks(dff,'MinPeakWidth',100,'MinPeakProminence',1.0*std(dff));
% 	hold on; plot(locs/framerate,pks+shiftup,'*','Color',c(i,:));% line([locs(:,:)/framerate,locs(:,:)/framerate],[-1 9],'Color',c(i,:))
% 	hold off
%     
%   locs = transpose(locs);
%    
%   hold on
%   line([locs(:,:)/framerate,locs(:,:)/framerate],[500 2000],'Color',c(i,:))
%   hold off