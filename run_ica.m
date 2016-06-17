%% LOAD VIDEO

% loading
load('../Recordings/SanneSleep23_LNY19RB_12-6-15_II_spont(20).mat');

%% CALCULATE ICA

% use nice wrapper to handle resizing and unrolling video
[ica_sig, mixing, separating, height, width] = NP_PerformICA(video.frames(30:end-30));

%% APPLY TO ORIGINAL VIDEO

% plot
plot_ica(ica_sig, mixing, height, width);

%% APPLY TO NEW VIDEO

% loading
load('../Recordings/SanneSleep23_LNY19RB_12-6-15_II_spont(22).mat');

% the actual ICA part:
% use separating matrix to perform saming unmixing
ica_sig2 = NP_ApplyICA(video.frames(30:end-30), separating, height, width);

% plot
plot_ica(ica_sig2, mixing, height, width);
