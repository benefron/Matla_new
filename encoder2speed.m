function [runningSpeed] = encoder2speed(data,r)
%A function to convert data from analog data of encoder into velocity
data(1:120000) = median(data);
data_rescale = rescale(data,-pi,pi);
data_unwrap = unwrap(data_rescale);
data_distance = r*data_unwrap;
figure; plot([1:18000000]/30000,data_distance(1:18000000));
[b,a] = butter(5,80/15000,'low');
data_filtered = filtfilt(b,a,data_distance);
runningSpeed = diff(data_filtered)*30000*100;
figure; plot([1:18000000]/30000,runningSpeed(1:18000000));
end

