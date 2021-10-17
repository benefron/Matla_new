function [data_binned_norm,time_vector] = psth(raster_mat,binSize)

sum_all = sum(raster_mat);
bin_to_samples = binSize * 30;
number_of_bins = 1:length(raster_mat)/bin_to_samples;
ind_mat(:,1) = 1 + number_of_bins * bin_to_samples - 300;
ind_mat(:,2) = number_of_bins * bin_to_samples;
data_binned = zeros(length(ind_mat),1);
for bin = 1: length(ind_mat)
    data_binned(bin) = sum(sum_all(ind_mat(bin,1):ind_mat(bin,2)));
end
data_binned_norm = data_binned/size(raster_mat,1);
window_size = length(raster_mat)/30000;
BF = length(number_of_bins)/window_size;
time_vector = [1:window_size*BF] - (window_size*(BF/2));
time_vector = time_vector/BF;
% figure;
% plot(time_vector,data_binned_norm)
% xlabel('time (sec)');
% ylabel('normalized firing rate (Hz)');

end
