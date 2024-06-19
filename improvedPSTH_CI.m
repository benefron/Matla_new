function [data_binned_norm_smoothed, time_vector, lowerCI, upperCI] = improvedPSTH_CI(raster_mat, binSize, varargin)
% Improved PSTH function with variable bin size and optional smoothing
% raster_mat: Matrix of raster data (neurons x time)
% binSize: Bin size in milliseconds
% varargin: Optional smoothing window size in bins

% Convert bin size from ms to samples (assuming 30 samples per ms)
bin_to_samples = round(binSize * 30);

% Calculate number of bins
number_of_bins = floor(size(raster_mat, 2) / bin_to_samples);

% Initialize index matrix for bins
ind_mat = zeros(number_of_bins, 2);
for i = 1:number_of_bins
    ind_mat(i, 1) = 1 + (i-1) * bin_to_samples;
    ind_mat(i, 2) = i * bin_to_samples;
end

% Calculate sum of all spikes in each bin for each trial
num_trials = size(raster_mat, 1);
spike_counts_per_bin = zeros(num_trials, number_of_bins);

for trial = 1:num_trials
    for bin = 1:number_of_bins
        spike_counts_per_bin(trial, bin) = sum(raster_mat(trial, ind_mat(bin, 1):ind_mat(bin, 2)));
    end
end

% Normalize binned data to get firing rate (Hz)
data_binned_norm = mean(spike_counts_per_bin, 1) / (binSize / 1000);

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

% Calculate 95% confidence interval
std_dev = std(spike_counts_per_bin / (binSize / 1000), 0, 1);
sem = std_dev / sqrt(num_trials); % Standard error of the mean
lowerCI = data_binned_norm - 1.96 * sem;
upperCI = data_binned_norm + 1.96 * sem;

% Optional: Plotting
% Uncomment the following lines to enable plotting within the function
% figure;
% plot(time_vector, data_binned_norm_smoothed);
% hold on;
% plot(time_vector, lowerCI, 'r--');
% plot(time_vector, upperCI, 'r--');
% xlabel('Time (s)');
% ylabel('Firing rate (Hz)');
% hold off;

end
