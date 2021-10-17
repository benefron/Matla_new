
function [Avisoft_times,aud_op] = avisoftSync(params_file)
%% load the files and create the variables
run(params_file);
data = load_open_ephys_binary([exp_path,'/experiment1/recording1/structure.oebin'],'continuous',1,'mmap');
channel_number = data.Header.num_channels;
aud_op = data.Data.Data.mapped(channel_number-(8-mic),:);
clear data
cd(folder_path);
files = dir('*.wav');
[aud_avi] = audioread([files.folder,'/',files.name]);

%% Find start point in both variables
aud_op = normalize(double(aud_op),"range",[-1,1]);
aud_avi = normalize(aud_avi,"range",[-1,1]);

avi_thr = iqr(aud_avi) * 5;
op_thr =  iqr(aud_op) * 5;

avi_cross = find(aud_avi > avi_thr,1) / 250000;
op_cross = find(aud_op > op_thr,1) / 30000;

% find differance in time before beep and shorten avisoft recordign

time_differance = floor((avi_cross - op_cross) * 250000);
aud_avi = aud_avi(time_differance:end);

% imprt artifact times
cd(exp_path)
load('condition_extract.mat');
artifacts_op = double(condition_extract.all_changes)-2*30000;
artifacts_op(1:end,2) = [artifacts_op(2:end,1);length(aud_op)];
% find artifact times for avisoft
artifacts_avi = floor(artifacts_op/30000*250000);

% Divide recordings based on artifact timings
data_cut = cell(2,length(artifacts_avi)+1);
data_cut{1,1} = aud_op(1:artifacts_op(1,1));
data_cut{2,1} = aud_avi(1:artifacts_avi(1,1));
for artifact = 1:length(artifacts_avi)
    data_cut{1,artifact+1} = aud_op(artifacts_op(artifact,1):artifacts_op(artifact,2));
    data_cut{2,artifact+1} = lag_corr(aud_avi(artifacts_avi(artifact,1):artifacts_avi(artifact,2)),data_cut{1,artifact+1});
end
clear aud_avi aud_op


[shortened_seg] = lag_corr(data_cut{2,10},data_cut{1,10});

%length_open = op_cross(end) - op_cross(1);
%time_total = length_open/SF;
%true_avi_SF = (avi_cross(end-1) - avi_cross(1))/time_total;

%SF_ratio = true_avi_SF/SF;
%time_vector_avi = [1:length(aud_avi)]' - avi_cross(1);
%time_vector_avi = time_vector_avi + (floor(SF_ratio * op_cross(1)));







%% create output structures
%Avisoft_times.avi_start = avi_cross(1);
%Avisoft_times.avi_end = avi_cross(end) + 1;
%Avisoft_times.alignedFrames = time_vector_avi;
%Avisoft_times.SF = true_avi_SF;
%Avisoft_times.open_ephys_start = op_cross(1);
%Avisoft_times.open_ephys_end = op_cross(end) + 1;
%Avisoft_times.sampling_ratio = SF_ratio;
end

