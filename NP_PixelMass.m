function [img_comp, mass_norm, std_proj_norm] = NP_PixelMass(video, varargin)
%NP_PIXELMASS Calculates the center of mass for each pixel over time
%   This function takes a (raw) video, performs basic background
%   subtraction on it, then calculates the center of mass for each pixel
%   over time. This was adapted from Will's FS_plot_allpxs function. The
%   function returns several variables for showing related data:
%
%   If you call the function without any outputs, then it will display the 
%   results in a new figure window. 
%
%   The first output, `img_comp`, is a RGB image showing the center of mass
%   for those pixels with a large std deviation.
%
%   The second output, `mass_norm` is the normalized center of mass for 
%   each pixel. 
%
%   The third output, `std_proj_norm` is normalized standard deviation for 
%   each pixel, used to highlight those pixels with activity of note. 
%
%   (The `img_comp` is a composite of the `mass_norm`
%   colored based on the jet colormap and the `std_proj_norm` as an alpha
%   channel.)

% parameters
filt_rad = 20; % gauss filter radius
filt_alpha = 30; % gauss filter alpha
lims = 3; % contrast prctile limits (i.e. clipping limits lims 1-lims)
cmap = jet(64);
per = 0; % baseline percentile (0 for min)

% load custom parameters
nparams=length(varargin);
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

% check type
if ~isa(video, 'double')
    video = mat2gray(video);
end

% get size
[sz_rows, sz_columns, sz_frames] = size(video);

% filter
h = fspecial('gaussian', filt_rad, filt_alpha);
video = imfilter(video, h, 'replicate');

% calculate baseline using percentile and repeat over video
baseline = prctile(video, per, 3);

% calculate dff
video_dff = bsxfun(@minus, video .^ 2, baseline .^ 2);
video_dff = bsxfun(@rdivide, video_dff, baseline);

% figure out center of max
com_idx = reshape(1:sz_frames, 1, 1, sz_frames);
com_idx = repmat(com_idx, sz_rows, sz_columns, 1);

% calculate mass
mass = sum(video_dff, 3);
mass = sum(video_dff .* com_idx, 3) ./ mass;

% normalize
mass_norm = mat2gray(mass, [1 sz_frames]);

% convert to color indices
idx_img = round(mass_norm .* size(cmap, 1));
img = ind2rgb(idx_img, cmap);

% calculate std dev projection
std_proj = std(video_dff, [], 3);

% normalize std projection
clims = prctile(std_proj(:), [lims 100 - lims]);
std_proj_norm = mat2gray(std_proj, clims);

% make composite image
img_comp = img .* repmat(std_proj_norm, 1, 1, 3);

if nargout == 0
    figure;
    imshow(img_comp);
    colormap(cmap); colorbar;
end

end

