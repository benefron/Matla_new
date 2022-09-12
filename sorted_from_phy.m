function [sorted_spikes] = sorted_from_phy(sorting_path)
%this function outputs a structure with 'good' and 'mua' fields containing
%the info about the sorted spikes
%   This function recives a folder path with the sorted data from kilosort2 and phy2
%   it reads the relevant files from the folder to create a structure
%   containing the data of spike timing divided into the appropriate units

cd(sorting_path)


%% read the relevant files into matlab format

% Read the .npy file with the spike times (in frame) of all units
spike_times = readNPY('spike_times.npy');

% Read the .npy file containing the cluster classification for each spike
spike_clusters = readNPY('spike_clusters.npy');

% Read a .tsv file containing additional information about the spikes and
% clusters (This code used the id field and the classification as 'good' or
% 'mua'
cluster_info = tdfread('cluster_info.tsv');


%% Identify the quality of each unit and allocate the spike times into the clusters

% Get the cluster ids for 'good' and 'mua'
try
    mua_cluster_id = cluster_info.cluster_id(all(ismember(cluster_info.group,'mua  '),2));
    good_cluster_id = cluster_info.cluster_id(all(ismember(cluster_info.group,'good '),2));
catch
    mua_cluster_id = cluster_info.id(all(ismember(cluster_info.group,'mua  '),2));
    good_cluster_id = cluster_info.id(all(ismember(cluster_info.group,'good '),2));
end
% creates an empty matrix to hold the spike times for each cluster
good_mat = zeros(length(good_cluster_id),max(cluster_info.n_spikes(all(ismember(cluster_info.group,'good '),2))));
mua_mat = zeros(length(mua_cluster_id),max(cluster_info.n_spikes(all(ismember(cluster_info.group,'mua  '),2))));

% Fill the empty 'good' mat with  the spike times of each cluster
for good = 1:length(good_cluster_id)
    temp_times = (spike_times(spike_clusters == good_cluster_id(good)))';
    good_mat(good,1:length(temp_times)) = temp_times;
end

% Fill the empty 'mua' mat with times of each cluster
for mua = 1:length(mua_cluster_id)
    temp_times = (spike_times(spike_clusters == mua_cluster_id(mua)))';
    mua_mat(mua,1:length(temp_times)) = temp_times;
end

% Creates the structure output of the function
sorted_spikes = struct;
sorted_spikes.good.id = good_cluster_id;
sorted_spikes.good.times = good_mat;
sorted_spikes.mua.id = mua_cluster_id;
sorted_spikes.mua.times = mua_mat;


