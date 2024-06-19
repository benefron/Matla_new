function [data_binned_norm_smoothed, time_vector] = improvedPSTH(raster_mat, binSize, varargin)
% Improved PSTH function with variable bin size and optional smoothing
% raster_mat: Matrix of raster data (neurons x time)
% binSize: Bin size in milliseconds
% varargin: Optional smoothing window size in bins

% Convert bin size from ms to samples (assuming 30 samples per ms)
bin_to_samples = round(binSize * 30);

% Calculate number of bins
number_of_bins = floor(size(raster_mat,2) / bin_to_samples);

% Initialize index matrix for bins
ind_mat = zeros(number_of_bins, 2);
for i = 1:number_of_bins
    ind_mat(i,1) = 1 + (i-1) * bin_to_samples;
    ind_mat(i,2) = i * bin_to_samples;
end

% Calculate sum of all spikes in each bin
sum_all = sum(raster_mat, 1);
data_binned = zeros(number_of_bins, 1);
for bin = 1:number_of_bins
    data_binned(bin) = sum(sum_all(ind_mat(bin,1):ind_mat(bin,2)));
end

% Normalize binned data to get firing rate (Hz)
data_binned_norm = data_binned / size(raster_mat,1) / (binSize / 1000);

% Check for smoothing argument
if ~isempty(varargin)
    window_size = varargin{1};
    % Smooth the PSTH using a moving average or another method
    % Example: Simple moving average smoothing
    data_binned_norm_smoothed = movmean(data_binned_norm, window_size);
else
    data_binned_norm_smoothed = data_binned_norm;
end

% Calculate time vector for plotting
time_vector = ((1:number_of_bins) - 0.5) * (binSize / 1000); % Center of each bin

% Optional: Plotting
% Uncomment the following lines to enable plotting within the function
% figure;
% plot(time_vector, data_binned_norm_smoothed);
% xlabel('Time (s)');
% ylabel('Firing rate (Hz)');

end
