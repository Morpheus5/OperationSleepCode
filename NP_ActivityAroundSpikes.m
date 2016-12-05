function NP_ActivityAroundSpikes(ROI, directory_raw, directory_out, varargin)
%NP_ACTIVITYAROUNDSPIKES
%   This function is meant to help with finding an extracting a window of
%   activity around spikes for comparison in a paper. This is achieved through
%   a two pass process. First, videos are analyzed once to figure out traces 
%   patterns (using NP_ExtractROIs) and spikes (using NP_From_Raw_To_Traces).
%   This is both done with donuts and without donuts. Those ROIs with a 
%   sufficient number of donut-spikes (specified by `min_donut_spikes`, which 
%   defaults to 7) are then used for the second stage.
%
%   During the second stage, the system will reload the videos and extract
%   frames in a window of time surrounding each spike for each ROI selected
%   during the first stage. A center of mass image, as well as a mat file,
%   are saved to the output directory (organized into subfolders for each
%   ROI).
%
%   There are a number of additional settings in the parameter section below
%   for customizing the behavior.

%% parameters
% video information
fs = 30;
skip_frames = 30; % first and last frames from each video that are discarded

% selection details
spike_donut_threshold = 3; % standard deviation for detecting spikes
spike_all_threshold = 3; % standard deviations for detecting spikes
min_donut_spikes = 7; % only generate figures for ROIs that have at least this many donut spikes

% what window to consider
window = [0.5 0.5]; % seconds before and after

% what to save
save_mat = true; % save a MATLAB file with traces, spikes, and other info
save_image = true; % save an un-annotated mass image
save_image_with_annot = true; % save a mass image with a circle around the ROI
save_plot = true; % save a plot of all traces
save_video = true; % save the video for the time slice

% other details
show_progress = true; % show a little progress bar
mask = '*.mat'; % mask for finding movies in the directory
debug = false; % open a window showing spike detection and other debugging

% load custom parameters
nparams = length(varargin);
if 0 < mod(nparams, 2)
	error('Parameters must be specified as parameter/value pairs');
end
for i = 1:2:nparams
    nm = lower(varargin{i});
    if ~exist(nm, 'var')
        error('Invalid parameter: %s.', nm);
    end
    eval([nm ' = varargin{i+1};']);
end

%% setup
% get list of files
files = dir(fullfile(directory_raw, mask));
files = cellfun(@(x) fullfile(directory_raw, x), {files(:).name}, 'UniformOutput', false);
files_count = length(files);

% convert window to samples
if isscalar(window)
    window = window * [1 1];
end
window = ceil(window * fs);

% convert skip frames to beginning / end
if isscalar(skip_frames)
    skip_frames = skip_frames * [1 1];
end

% open progress bar
if show_progress
    h = waitbar(0, 'Extracting data...');
end

%% find ROIs of interest
% this is the first pass through the videos; loading the videos twice is
% pretty slow, but it maximizes the

spikes_donut = cell(size(files)); % ROI * time
spikes_all = cell(size(files)); % ROI * time
traces_all = cell(size(files)); % ROI * time

for i = 1:files_count
    % load video (slow)
    v = load(files{i}, 'video*');
    nms = fieldnames(v);
    if 1 ~= length(nms)
        warning('More than one variable in the video file. Skipping.');
        continue;
    end
    if isfield(v.(nms{1}), 'frames')
        video = v.(nms{1});
    else
        video = struct('frames', v.(nms{1}));
    end
    clear v;
    
    % convert to grayscale
    if 3 == ndims(video.frames(1).cdata) && 3 == size(video.frames(1).cdata, 3)
        video_gs = video_rgb2gray(cat(4, video.frames(:).cdata));
        video.frames = struct('cdata', squeeze(num2cell(video_gs, [1 2])));
        clear video_gs;
    end
    
    % extract ROIs with donut
    roi_trace = NP_ExtractROIs(video.frames((1 + skip_frames(1)):(end - skip_frames(2))), ROI, 1);
    
    % convert to spikes
    [~, ~, roi_spikes] = NP_From_Raw_To_Traces({roi_trace}, spike_donut_threshold);
    
    % store spikes
    spikes_donut{i} = squeeze(roi_spikes);
    
    % extract ROIs without donut
    roi_trace = NP_ExtractROIs(video.frames((1 + skip_frames(1)):(end - skip_frames(2))), ROI, 0);
    
    % convert to spikes
    [~, ~, roi_spikes] = NP_From_Raw_To_Traces({roi_trace}, spike_all_threshold);
    
    % store spikes
    traces_all{i} = roi_trace;
    spikes_all{i} = squeeze(roi_spikes);
    
    % clear video
    clear video;
    
    % if debug
    if debug
        figure;
        subplot(1, 3, 1);
        plot_many(traces_all{i}');
        title(files{i});
        subplot(1, 3, 2);
        plot_many(spikes_all{i}');
        title('All spikes');
        subplot(1, 3, 3);
        plot_many(spikes_donut{i}');
        title('Donut spikes');
    end
    
    % update progress
    if show_progress
        waitbar(i / (2 * files_count));
    end
end

%% identify ROIs of interest
donut_spikes_per_roi = sum(cat(2, spikes_donut{:}), 2);

% rois of interest
rois_of_interest = donut_spikes_per_roi >= min_donut_spikes;
rois_of_interest = reshape(rois_of_interest, 1, []); % make row, better for iteration

% make folders for each ROI
for i = find(rois_of_interest)
    fprintf('ROI %d: %d donut spikes\n', i, donut_spikes_per_roi(i));
    if ~exist(fullfile(directory_out, num2str(i)), 'dir')
        mkdir(directory_out, num2str(i));
    end
end

%% using ROIs of interest, look through videos to find surrounding vehavior

% reusable mesh for drawing rois
x = []; y = [];
video_width = 0; video_height = 0;

for i = 1:files_count
    % get spikes and traces for current video
    cur_spikes = spikes_all{i};
    cur_traces = traces_all{i};
    
    % skipped file
    if isempty(cur_spikes)
        continue;
    end
    
    % has any spikes of interest? if not, skip
    if ~any(any(cur_spikes(rois_of_interest, (1 + window(1)):(end - window(2)))))
        continue;
    end
    
    % load video (slow)
    v = load(files{i}, 'video*');
    nms = fieldnames(v);
    if 1 ~= length(nms)
        warning('More than one variable in the video file. Skipping.');
        continue;
    end
    if isfield(v.(nms{1}), 'frames')
        video = v.(nms{1});
    else
        video = struct('frames', v.(nms{1}));
    end
    clear v;
    
    % convert to grayscale
    if 3 == ndims(video.frames(1).cdata) && 3 == size(video.frames(1).cdata, 3)
        video_gs = video_rgb2gray(cat(4, video.frames(:).cdata));
        video.frames = struct('cdata', squeeze(num2cell(video_gs, [1 2])));
        clear video_gs;
    end
    
    % name
    [~, name, ~] = fileparts(files{i});
    
    % for each roi
    for j = find(rois_of_interest)
        % for each spike
        for k = find(cur_spikes(j, :))
            % make sure window is available
            if k < (1 + window(1))
                continue;
            end
            if k > (size(cur_spikes, 2) - window(2))
                continue;
            end
            
            % get relevant frames
            frames = cat(3, video.frames((skip_frames(1) + k - window(1)):(skip_frames(1) + k + window(2))).cdata);
            % scaling greyscale 
            raw = video_adjust(frames, [0.5 0.999]);
            
            % extract traces
            traces = cur_traces(:, (k - window(1)):(k + window(2)));
            spikes = cur_spikes(:, (k - window(1)):(k + window(2)));
            
            % save plot
            if save_plot
                f = figure;
                plot_many(((skip_frames(1) + k - window(1)):(skip_frames(1) + k + window(2))) ./ fs, traces');
                title(sprintf('%s; ROI %d', strrep(name, '_', '\_'), j));
                print(gcf, fullfile(directory_out, num2str(j), sprintf('%s_%d.png', name, k)), '-dpng', '-r300');
                close(f);
            end
                        
            % calculate mass image
            cmap = jet(64);
            % use this line to compress the colormap to show variation that
            % is close in time (on a small spectrum of the colormap):
            cmap = [repmat(cmap(1, :), 16, 1); cmap; repmat(cmap(end, :), 16, 1)];
            % then run NP_PixelMass
            [im, im_mass, im_std] = NP_PixelMass(frames, 'cmap', cmap);
            if save_image
                imwrite(im, fullfile(directory_out, num2str(j), sprintf('%s_%d.jpg', name, k)));
            end

            if save_image_with_annot
                % check mesh grid
                if size(frames, 1) ~= video_height || size(frames, 2) ~= video_width
                    video_height = size(frames, 1);
                    video_width = size(frames, 2);
                    [x, y] = meshgrid(1:video_width, 1:video_height);
                end
                
                % stats
                center = ROI.stats(j).Centroid;
                radius2 = (ROI.stats(j).Diameter / 2) ^ 2;

                % calculate distance
                distance = (x - center(1)) .^ 2 + (y - center(2)) .^ 2;
                mask_space = distance > radius2 & distance < radius2 * 2;
                im(repmat(mask_space, 1, 1, 3)) = 1;
                imwrite(im, fullfile(directory_out, num2str(j), sprintf('%s_%d_self.jpg', name, k)));
            end
            
            if save_video
                % dff
                dff = NP_Dff(frames);
                
                % draw spikes
                d = false(size(spikes));
                d(j, :) = spikes(j, :);
                dff2 = NP_DrawROIs(dff, ROI, d, fs);
                
                % write
                video_write(fullfile(directory_out, num2str(j), sprintf('%s_%d_raw.mp4', name, k)), raw, fs); % black screen for sleep videos
                %video_write(fullfile(directory_out, num2str(j), sprintf('%s_%d_dff.mp4', name, k)), dff, fs);
                %video_write(fullfile(directory_out, num2str(j), sprintf('%s_%d_WithSpikes.mp4', name, k)), dff2, fs);
            end
            
            if save_mat            
                % meta data for saving
                roi = j;
                file = files{i};
                frame = skip_frames(1) + k;

                % save matlab data
                save(fullfile(directory_out, num2str(j), sprintf('%s_%d.mat', name, k)), ...
                    'roi', 'file', 'frame', 'traces', 'spikes', 'im_mass', 'im_std');
            end
        end
    end
    
    % clear video
    clear video;
    
    % update progress
    if show_progress
        waitbar((files_count + i) / (2 * files_count));
    end
end

%% clean up
% close progress
if show_progress
    close(h);
end

end
