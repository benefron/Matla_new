function [times] = event_end_detection(data,timeBefore)


%sigmaTimeSec = 0.08;
%sigmaSamples = round(sigmaTimeSec*30000);
%window = fspecial('gaussian',[1, sigmaSamples*6],sigmaSamples);
%data = conv(data,window,'same');
data = abs(data);

figure;
plot(data);
yline(iqr(data))
thr = input(['black line = ' num2str(iqr(data)),'; input threshold value: '])
close all
figure;
plot(data)
yline(thr);
thr = input('reenter threshold or correct')


above_thr = find(data > thr)+1;
below_thr = find(data < thr);
times = intersect(above_thr',below_thr');


times((times - timeBefore < 1)) = [];
try
    wind_mat = (ones(length(times),timeBefore) .* [1:1:timeBefore]) + times;
catch
    wind_mat = (ones(length(times),timeBefore) .* [1:1:timeBefore]) + times';
end
wind_mat = round(wind_mat);
wind_dat_onset = data(wind_mat);
[rows,~] = find(wind_dat_onset < (thr)); rows = unique(rows);
times(rows) = [];

