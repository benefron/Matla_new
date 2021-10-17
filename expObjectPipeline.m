%% Create a auditory ephys experiments object and run analysis


%% sort spikes


%% create the experiment object

parameters_file_path = []; %Input the parameters file path that contain the experiment data

experiement = auditory_ephys_exp(parameters_file_path); % substitute experiment with the appropriate experiment name 

%% Sync the avisoft data to the ephys data

 % Based on mlx file in FVB/757/audio
 
 % analyze audio and extract aluminum events using avisoft lab software
 
 % add the start and end times of the aluminum events to the main object
 