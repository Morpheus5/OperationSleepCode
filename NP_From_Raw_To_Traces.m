function [Av_Data,Trace_Data,Spike_Data] = NP_From_Raw_To_Traces(aveRoiRaw,threshold,Fs,binSize)
% Created: 2015/12/11 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/12/11
% By: Byron Price
%   This file will run through raw data of fluorescence traces and output
%    detrended fluorescence traces, spike data, and summed spike data
%    across ROIs (called avalanches).
% Updated: 2016/05/01
% By: Nathan Perkins
%   Eliminated non-linear fitting, consolidated spike detection,
%   introduced more flexibility in spike detection criteria (including
%   requiring consecutive moments above threshold) and support finding when
%   the spike begins.
%
% INPUT:   aveRoiRaw - the roi_ave.raw matrix from a single session
%          FOUR OPTIONAL INPUTS
%          threshold - z-score threshold (peaks in activity above this value
%                       are counted as spikes, default = 2
%          Fs - sampling frequency (Hz), default = 30
%          binSize - width of bins (secs), spikes are counted as having occurred
%                       within bins of a certain size, default = 1/Fs
% OUTPUT:  Av_Data - number videos by number frames matrix of summed spike
%                       activity
%          Trace_Data - number videos by number ROIs by number frames
%                       matrix of detrended fluorescence traces
%          Spike_Data - number videos by number ROIs by number frames
%                       matrix of spike data (zeros and ones for each ROI
%                       indicating a spike in a bin)
%

tic

%OPTIONAL INPUTS
if nargin > 5
    error('Gardner:DataStream:TooManyInputs', ...
        'Requires at most 4 optional inputs');
end

% Fill in unset optional inputs.
switch nargin
    case 1
        threshold = 3;
        Fs = 30;
        binSize = 1/Fs;
    case 2
        Fs = 30;
        binSize = 1/Fs;
    case 3
        binSize = 1/Fs;
end

% FIGURE OUT HOW MANY VIDEOS, ROIs, and FRAMES
alldata = aveRoiRaw;
numVideos = size(alldata,2);
numROIs = size(alldata{1},1);
numFrames = length(alldata{1,1});
consecutive = 2;
find_start = true;

% DETREND EACH ROI INDIVIDUALLY FOR EACH VIDEO
Trace_Data = zeros(numVideos,numROIs,numFrames);
for i = 1:numVideos
    for j = 1:numROIs
        [alldata{1,i}(j,:)] = Preprocessing(alldata{1,i}(j,:));
    end
    Trace_Data(i,:,:) = alldata{1,i};
end

Detrended = toc

% COUNT SPIKES AND SUM ACROSS ROIs IN A SINGLE RECORDING
Spike_Data = false(numVideos,numROIs,numFrames);
Av_Data = zeros(numVideos,numFrames);
for j=1:numVideos
    binarySpikes = []; % will contain event data for all ROIs in 
                        % a single recording
    for k=1:numROIs
        % detect calcium spikes
        [~,Spikes] = Spike_Detector(squeeze(Trace_Data(j,k,:)),Fs,binSize,threshold,consecutive,find_start);
        binarySpikes = [binarySpikes,Spikes];
    end
    Spike_Data(j,:,:) = binarySpikes';
    sumSpikes = sum(binarySpikes,2); % add up spikes from each ROI
                                        % into a single vector, this
                                        % will be used to deterimine
                                        % avalanche size

    Av_Data(j,:) = sumSpikes;
           
end

save('Spontaneous_Data.mat','Av_Data','Spike_Data','Trace_Data');

end



function [timeSeries] = Preprocessing(timeSeries)
% Av_Preprocessing.m
%   Detrend a time series and subtract out the mean
% Created:2015/10/21 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/12/11
%   By: Byron Price
% INPUT: timeSeries - calcium imaging pixel intensities over time for a
%                       single ROI
% OUTPUT: timeSeries - detrended data

% LINEAR REGRESSION TO SUBTRACT OUT BLEACHING
x = 1:1:length(timeSeries);
X = [ones(length(x),1) x'];
b = X\timeSeries';
ycalc = X*b;

timeSeries = (timeSeries'-ycalc)';

timeSeries = timeSeries./max(timeSeries);
end

function [numBins,Spikes] = Spike_Detector(timeSeries,sampleFreq,binSize,threshold, consecutive, find_start)
%Avalanche.m  
%   Detect calcium spikes in imaging data for a single ROI.
% Created: 2015/09/30 at 24 Cummington, Boston
%  Byron Price
% Updated: 2015/12/11
%  By: Byron Price
% Updated: 2016/05/01
%  By: Nathan Perkins
%
% REFERENCE: Klaus, Plenz 2011 Statistical Analyses Support Power Law ...
%            Ribeiro, Copelli et al. 2010 Spike Avalanches Exhibit ...
%
% INPUT: timeSeries - change in activity over time for a single ROI in
%               units given by T = 1/Fs  ... Fs = sampling frequency
%        sampleFreq - sampling frequency, Hz
%        binSize - width of bin (seconds)
%               Spikes are counted as having occurred within a given 
%               bin if the activity of that bin is above "threshold"
%        threshold - zscore threshold (1.5, 2 etc. standard deviations
%               above the mean) for decision: above threshold = activity 
%               forms part of an avalanche
%        consecutive - number of consecutive timesteps above threshold
%                   required to count as spike
%        find_start - whether or not to find the start
% OUTPUT: numBins - given the binSize and the length of the recording,
%               numBins = (total # frames) / (# frames / bin)
%         Spikes - a vector containing either a 0 or a 1, 1 being an
%               spike of activity within that bin, 0 otherwise

T = 1/sampleFreq; % period for a single frame
framesPerBin = round(binSize/T); 


times = 1:framesPerBin:length(timeSeries);
numBins = length(times); %number of bins, based on the
                         % length of the data and bin size (secs)

% timeSeries = timeSeries./max(timeSeries);
stdActivity = std(timeSeries);

% SPIKE DETECTION

Spikes = false(size(timeSeries));

% the trace was already divided by its maximum value, there are a few steps to being
% counted as a calcium spike:

% 1) have value greater than `threshold` standard deviations above mean
% (defaults to 2 standard deviations)
above_threshold = (timeSeries ./ stdActivity) > threshold;

% calculate deltas
timeSeriesDelta = [0; diff(timeSeries)];

% find starts of spikes
start_above_threshold = above_threshold & [true; ~above_threshold(1:(end-1))];
for j = find(start_above_threshold')
    % find length of time above threshold
    l = find(~above_threshold(j:end), 1);
    if isempty(l)
        % no moment below threshold? remainder of spike train is above threshold
        l = 1 + length(above_threshold) - j;
    else
        l = l - 1;
    end
    
    % 2) require more than `consecutive` moments above the threshold
    if l < consecutive
        continue;
    end
    
    % 3) find the local maximum
    [~, spike_peak] = max(timeSeries(j:(j + l - 1)));
    spike_peak = spike_peak + j - 1;
    
    % 4) step back to "start" of spike
    if find_start
        spike_start = spike_peak;
        while spike_start > 1 && timeSeriesDelta(spike_start - 1) >= 0
            spike_start = spike_start - 1;
        end
        Spikes(spike_start) = true;
    else
        Spikes(spike_peak) = true;
    end
end

end  


