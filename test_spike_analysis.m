%% synthetic data
roi_count = length(ROI.stats);
sample_count = 2000;

% STIMULUS: random
% Spike_Data = rand(roi_count, sample_count) < 0.01;
% Spike_Data = reshape(Spike_Data, 1, roi_count, sample_count);

% STIMULUS: sequential, non-overlapping
% space = 0;
% Spike_Data = false(roi_count, sample_count);
% for i = 1:roi_count
%     st = 1 + (i - 1) * (1 + space);
%     ev = (1 + space) * roi_count;
%     Spike_Data(i, st:ev:end) = true;
% end
% Spike_Data = reshape(Spike_Data, 1, roi_count, sample_count);

% STIMULUS: sequential, overlapping
% space = 4;
% Spike_Data = false(roi_count, sample_count);
% for i = 1:roi_count
%     st = i;
%     ev = (1 + space);
%     Spike_Data(i, st:ev:end) = true;
% end
% Spike_Data = reshape(Spike_Data, 1, roi_count, sample_count);

% STIMULUS: actual
% load('../Recordings/roi/Spontaneous_Data.mat');

% STIMULUS: all co-active
% Spike_Data = false(roi_count, sample_count);
% every = 20;
% Spike_Data(:, every:every:end) = true;
% Spike_Data = reshape(Spike_Data, 1, roi_count, sample_count);

% STIMULUS: artificial 1 leads all others
every = 50;
Spike_Data = false(roi_count, sample_count);
Spike_Data(1, 1:every:end) = true;
Spike_Data(2:2:roi_count, 5:every:end) = true;
Spike_Data = reshape(Spike_Data, 1, roi_count, sample_count);

%% run
SM_SortedSpikeAnalysis

%% plot spikes
figure;
imagesc(squeeze(Spike_Data(1, :, 1:200)));
