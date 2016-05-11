%% BACKGROUND SUBTRACTION, NO BAGELS

% make dff
video_dff = NP_Dff(video.frames(30:end-30), 'clim', [0 1]);

% extract ROIs
ave_roi = NP_ExtractROIs(video_dff, ROI, 0); % 0 means no bagsels or donuts

% extract spikes
[~, ~, spikes] = NP_From_Raw_To_Traces({ave_roi}, 2);

%% NO BACKGROUND SUBTRACTION, YES BAGELS

% extract ROI
ave_roi = NP_ExtractROIs(video.frames(30:end-30), ROI, 1); % 1 means yes bagels

% extract spikes
[~, ~, spikes] = NP_From_Raw_To_Traces({ave_roi}, 2);

%% NO BACKGROUND SUBTRACTION, YES BAGELS -- "full version"

close all; clear all
load('LNY20RB_RH_9-7-15_SanneSleep6(6).mat')
load('image_roi/roi_data_image.mat');
ave_roi = NP_ExtractROIs(video.frames(30:end-30), ROI, 1); % 1 means yes bagels
[Av_Data, Trace_Data, Spike_Data] = NP_From_Raw_To_Traces({ave_roi}, 2);
SM_SortedSpikeAnalysis
save('Video6_ave_roi_SpikeData_NoBkgrSubYesBagels')