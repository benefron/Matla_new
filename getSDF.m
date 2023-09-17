function spikeVec = getSDF(unitTrain,sigma)

%sigma = 150;
edges = [-3*sigma:1:3*sigma];
kernel = normpdf(edges,0,sigma);
s = conv(unitTrain,kernel);
center = ceil(length(edges)/2);
s = s(center:end-center-1);
spikeVec = resample(s,1000,30000);
end