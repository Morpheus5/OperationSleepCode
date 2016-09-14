% load in the normalized center of mass per pixel matrices (512x512 double)  
% saved from NP_ActivityAroundSpikes 

% % I used sumv as control (similarity across cells as baseline dissimilarity)

cellnr = 1; % adjust so that the right traces are selected for plotting
DIR = pwd;
sumvlist = dir(fullfile(DIR,'*.mat')); % use for sleeping data
sumvlist = {sumvlist(:).name};
NrSpikes = length(sumvlist);

% prepare matrices for pixel mass dissimilarity scores:
total_array = [];
d = zeros(NrSpikes, NrSpikes);
sumd = [];

% prepare matrices for ROI trace dissimilarity scores:
total_traces = [];
d2 = zeros(NrSpikes, NrSpikes);
sumd2 = [];

% Combine all mass matrices into one array:
for iter = 1:NrSpikes
    [~,file,~] = fileparts(sumvlist{iter});
  	load(fullfile(DIR,sumvlist{iter}), 'im_*');
    load(fullfile(DIR,sumvlist{iter}), 'traces');
    
    % now that we know the image dimensions, setup our matrices
    if iter == 1
        total_array = zeros(NrSpikes, size(im_mass, 1), size(im_mass, 2));
        total_traces = zeros(NrSpikes, size(traces, 1), size(traces, 2));
    end
    
    % add image to array
    total_array(iter, :, :) = im_mass .* im_std;
    total_traces(iter, :, :) = traces;
end

% make difference matrix 'd' & summed dissimilarity vector 'sumd':
for i = 1:NrSpikes
   for j = 1:NrSpikes
       if j>i;
           difference = total_array(i, :, :) - total_array(j, :, :);
           d(i, j) = mean(abs(difference(:))); % MAE: mean absolute error
           sumd(end+1) = d(i, j);
           
           difference2 = total_traces(i, :, :) - total_traces(j, :, :);
           d2(i, j) = mean(abs(difference2(:))); 
           sumd2(end+1) = d2(i, j);
       end
    end
end

% save summed dissimilarity vector for all sleep cells:
save('Dissimilarity_matrix','sumd', 'sumd2','total_traces','-v7.3')

% plot summed dissimilarity distribution:
figure(1); hold on; title('Dissimilarity histogram pixel mass');
hist(sumd, 200); hold off
saveas(1,'Dissimilarity Histogram Pixel Mass.eps')

figure(2); hold on; title('Dissimilarity histogram traces');
hist(sumd2, 200); hold off
saveas(2,'Dissimilarity Histogram Traces.eps')

tracedata = squeeze(total_traces(:,cellnr,:));
figure(3); plot_many(tracedata')
saveas(3,'Spike Traces.eps')

% %% dissimilarity between-cells as control
% 
% % Make sure sumv singing and sumv sleep are loaded into the workspace
% 
% control_sumv_sleep = sumv_out_sleep;
% [a,b] = size(control_sumv_sleep);
%     
% for i = 1:a
%     tt = control_sumv_sleep(i,:);
%     tt(9900:11000) = mean(tt);
%     tt = tt-mean(tt);
%     control_sumv_sleep(i,:) = (tt-mean(tt))/(std(tt));
% end
% 
% d = zeros(a);
% for i = 1:a
%     for j = 1:a 
%         d(i,j) = sum(abs(control_sumv_sleep(i,:)-control_sumv_sleep(j,:)));
%     end
% end
% 
% % make summed dissimilarity vector:
% sumd_control_sleep = [d(1:length(d)*length(d))];
% 
% % save summed dissimilarity vector for all sleep cells:
% save('Control_Sleep_RA_dissimilarity_matrix','sumd_control_sleep')
% 
% % plot summed dissimilarity distribution:
% figure()
% hist(sumd_control_sleep,200)
% xlim([1 3e4]);
