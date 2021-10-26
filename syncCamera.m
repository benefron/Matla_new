function [cams] = syncCamera(exp_path,x,strobe_channel)
%Sync the camera frame times to the open ephys strobe times
%   Choose the camera file for face camera or for whisker camera and
%   extract the frames, find the starting frame from both the camera csv
%   file and the strobe data from the open ephys. 
%   x - categorical. recieves 1 or 2, 1 for whisker camera, 2 for face cam
% Choose if whisking camera or face camera and choose the appropriate file
cd(exp_path)
switch x
    case 1
        disp('Choose csv file for the whisking camera');
        %strobe_channel = Whiskers_cam;
    case 2
        disp('Choose csv file for the face camera');
        %strobe_channel = face_camera;
end
[file,path] = uigetfile('*.csv');
cam_fps = input(['Please input the fps for the current cam' '\n']);


% read csv file to variable
csv_data = csvread([path,'/',file]);
% read openephys all timestamps and digital channels timestamps npy files
full_timestamps = readNPY('experiment1/recording1/continuous/Rhythm_FPGA-153.0/timestamps.npy');
open_ephys_TTL = load_open_ephys_binary('experiment1/recording1/structure.oebin','events',1);
% find camera strobe on times from timestamps
strobe_data = double(open_ephys_TTL.Timestamps(open_ephys_TTL.Data == strobe_channel));
% find the frames from both openephys and csv after the forced stop
strobe_data(2:end,2) =  (strobe_data(2:end,1) - strobe_data(1:end-1))/open_ephys_TTL.Header.sample_rate;
csv_data(2:end,3) = (csv_data(2:end,2) - csv_data(1:end-1,2))/10000000;
strobe_start = find(strobe_data(1:3*cam_fps,2) > (5*(1/cam_fps)));
strobe_start = strobe_data(strobe_start(end),1);
csv_start = find(csv_data(1:3*cam_fps,3) > (5*(1/cam_fps)));
csv_start = csv_start(end);
% align both time stamps
csv_times = (csv_data(:,2) - csv_data(1,2))/10000000;
csv_times = csv_times - csv_times(csv_start);
csv_frames = floor(csv_times * open_ephys_TTL.Header.sample_rate);
csv_frames = csv_frames + (strobe_start - double(full_timestamps(1)));
% save the timestamps and frame times
cams.strobe_start_frame = strobe_start - double(full_timestamps(1));
cams.csv_start_frame = csv_start;
cams.csv_aligned_frames = csv_frames;
cams.fps = cam_fps;
end

