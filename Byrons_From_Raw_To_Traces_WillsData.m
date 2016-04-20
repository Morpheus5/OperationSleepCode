function [Av_Data,Trace_Data,Spike_Data] = Byrons_From_Raw_To_Traces_WillsData(AllWillsData,threshold,FitType,Fs,binSize)
%From_Raw_To_Traces.m 
%  
% Created: 2015/12/11 at 24 Cummington, Boston
%   Byron Price
% Updated: 2015/12/11
% By: Byron Price
%   This file will run through raw data of fluorescence traces and output
%    detrended fluorescence traces, spike data, and summed spike data
%    across ROIs (called avalanches).
%
% INPUT:   aveRoiRaw - the roi_ave.raw matrix from a single session
%          FOUR OPTIONAL INPUTS
%          FitType - string with type of fit to perform on data from single
%                       ROIs within individual recordings, this is for 
%                       detrending bleaching and some noise, MATLAB allows
%                       varied input, like 'poly1' for linear or 'fourier8'
%                       for 8-term Fourier series
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
        threshold = 2;
        FitType = 'fourier5'; % could try 'fourier1' 'fourier4' fourier8' 'exp1'
        Fs = 30;
        binSize = 1/Fs;
    case 2
        FitType = 'fourier5';
        Fs = 30;
        binSize = 1/Fs;
    case 3
        Fs = 30;
        binSize = 1/Fs;
    case 4
        binSize = 1/Fs;
end

% FIGURE OUT HOW MANY VIDEOS, ROIs, and FRAMES
data = AllWillsData; 
numVideos = size(data,2);
numROIs = size(data{1},2);
shift = 0; 
numFrames = size(data{1},1);

% DETREND EACH ROI INDIVIDUALLY FOR EACH VIDEO
Trace_Data = zeros(numVideos,numROIs,numFrames);
Data = zeros(numVideos,numFrames,numROIs);
Flipped_Data = zeros(numVideos,numROIs,numFrames);
alldata = cell(1,numVideos);

for ii = 1:numVideos
    Data(ii,:,:) = data{1,ii};
    Flipped_Data(ii,:,:) = permute(Data(ii,:,:),[1 3 2]);
    alldata{1,ii}(:,:) = Flipped_Data(ii,:,:);
end

for i = 1:numVideos
    for j = 1:numROIs
        [alldata{1,i}(j,:)] = Preprocessing(alldata{1,i}(j,:),FitType);
    end
    Trace_Data(i,:,:) = alldata{1,i}(:,:);
end

Detrended = toc

% COUNT SPIKES AND SUM ACROSS ROIs IN A SINGLE RECORDING
Spike_Data = zeros(numVideos,numROIs,numFrames);
Av_Data = zeros(numVideos,numFrames);
for j=1:numVideos
    binarySpikes = []; % will contain event data for all ROIs in 
                        % a single recording
    for k=1:numROIs
        % detect calcium spikes
        [~,Spikes] = Spike_Detector(squeeze(Trace_Data(j,k,:)),Fs,binSize,threshold);
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



function [timeSeries] = Preprocessing(timeSeries,FitType)
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
% x = 1:1:length(timeSeries);
% X = [ones(length(x),1) x'];
% b = X\timeSeries';
% ycalc = X*b;
% 
% timeSeries = (timeSeries'-ycalc)';

% REGRESSION TO SUBTRACT OUT BLEACHING AND OTHER 
% INHOMOGENEITIES IN THE DATA
x = 1:length(timeSeries);
f = fit(x',timeSeries',FitType);
Y = f(x);

timeSeries = (timeSeries'-Y)';

timeSeries = timeSeries./max(timeSeries);
end

function [numBins,Spikes] = Spike_Detector(timeSeries,sampleFreq,binSize,threshold)
%Avalanche.m  
%   Detect calcium spikes in imaging data for a single ROI.
% Created: 2015/09/30 at 24 Cummington, Boston
%  Byron Price
% Updated: 2015/12/11
%  By: Byron Price
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
% OUTPUT: numBins - given the binSize and the length of the recording,
%               numBins = (total # frames) / (# frames / bin)
%         Spikes - a vector containing either a 0 or a 1, 1 being an
%               spike of activity within that bin, 0 otherwise

T = 1/sampleFreq; % period for a single frame
framesPerBin = round(binSize/T); 


times = 1:framesPerBin:length(timeSeries);
numBins = length(times); %number of bins, based on the
                         % length of the data and bin size (secs)

Spikes = zeros(numBins,1); 
% timeSeries = timeSeries./max(timeSeries);
stdActivity = std(timeSeries);
% SPIKE DETECTION LOOP

% the trace was already divided by its maximum value, there are two steps to being
% counted as a calcium spike: 1) have value that is greater than or equal to
% 2 standard deviations above the mean; and 2) be a local maximum, meaning
% the fluorescence intensity of the current time step should be greater
% than both of its neighbors

for i=2:numBins-1
        if timeSeries(i)/stdActivity > threshold % tentatively counted as spike if above threshold
            if timeSeries(i) > timeSeries(i-1) && timeSeries(i) > timeSeries(i+1)
                Spikes(i) = 1;
            end
        end
end
end  


