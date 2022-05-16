function [times] = event_detection(data,timeBefore)


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


above_thr = find(data > thr);
below_thr = find(data < thr)+1;
times = intersect(above_thr',below_thr');


times((times - timeBefore < 1)) = [];
try
    wind_mat = (ones(length(times),timeBefore) .* [-timeBefore:1:-1]) + times;
catch
    wind_mat = (ones(length(times),timeBefore) .* [-timeBefore:1:-1]) + times';
end
wind_mat = round(wind_mat);
wind_dat_onset = data(wind_mat);
[rows,~] = find(wind_dat_onset > (thr)); rows = unique(rows);
times(rows) = [];



to_plot = 0; % If you want to varify the detection on a plot make 1
switch to_plot
    case 0
    case 1
        figure; plot ((data),'LineWidth',0.5);
        hold on
        for i=1:45%length(whisk_onset)
            xline(times(i),'r');
        end
        g = ones(length(times(5:30)),60000);
        g = g.*(-29999:30000);
        gx = g + times(5:30);
        gd = data(gx);
        figure; hold on; plot(mean(gd))
        %x = 1;
end




end

