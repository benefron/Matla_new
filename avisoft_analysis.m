function [cleaned_sound,aluminumCondition_sound] = avisoft_analysis(params)
run(params);
avi_times = avisoftSync(params_path);


cd(folder_path);
filename = dir('*.wav');

[data] = audioread(filename.name);
fs = avi_times.SF;

total_length = length(data);
seg_length = floor(total_length/20);
start = 0;  
finish = seg_length;
data_seg = zeros(seg_length,19);
data_last = data(seg_length*19+1:end);

for seg = 1:19
    data_seg(:,seg) = data(start+1:finish);
    start = finish;
    finish = finish + seg_length;
end
clear data
data_last = wdenoise(data_last,4);
parfor(dns = 1:19,2)
    data_seg(:,dns) = wdenoise(data_seg(:,dns),4);
end

cleaned_sound = reshape(data_seg,[],1);
cleaned_sound = [cleaned_sound;data_last];
clear data_last data_seg

cd(exp_path)
load('condition_extract.mat')

resampled_changes = condition_extract.all_changes*avi_times.sampling_ratio;
time_before = floor(condition_extract.artifact_times(1)/1000 * fs);
time_after = floor(condition_extract.artifact_times(2)/1000 * fs);
all_rmv = ones(time_after+time_before,length(condition_extract.all_changes));
rmv_vec = [-time_before:time_after-1]';
all_rmv = all_rmv .* rmv_vec;

all_rmv = all_rmv + double(resampled_changes)';

cleaned_sound(all_rmv) = 0;


