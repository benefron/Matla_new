function [spikeVec,unitBinned] = getSDFBinned(unitTrain,expLength)

% adapted from Matlab for Neuroscience
bin_edges = [0:0.001:expLength+0.001];
unitBinned = histcounts(unitTrain,bin_edges);
sigma = 0.015;
edges = [-3*sigma:0.001:3*sigma];
kernel = normpdf(edges,0,sigma);
kernel = kernel.*0.001;
s = conv(unitBinned,kernel);
center = ceil(length(edges)/2);
s = s(center:end-(center-1)); 
spikeVec = s;
% spikeVec = resample(s,1000,30000);
end