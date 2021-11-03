function [] = getConditions(ParametersFile)
run(ParametersFile); % runs .m file with the exp peramiteres, channel number for different inputs, pathways for the files etc.
cd(exp_path)
open_ephys_TTL = load_open_ephys_binary('experiment1/recording1/structure.oebin','events',1);
channel_states = open_ephys_TTL.Data;
timestamps = open_ephys_TTL.Timestamps;
cd([exp_path,'/experiment',num2str(recording_num),'/recording1/continuous/',FPGA_folder])
all_timestamps = readNPY('timestamps.npy');

timestamps = timestamps - all_timestamps(1);
OnTimestamps = timestamps(channel_states == arduino_channel);

isLastOn = [(OnTimestamps(2:end) - OnTimestamps(1:end-1)) > ((length_on_off * 2.5) * (SF/1000));1];
positions = find(isLastOn);
condition_classification = [positions(1);(positions(2:end) - positions(1:end-1))];
all_changes = OnTimestamps(positions);

condition_extract.all_changes = all_changes;
condition_extract.condition_classification = condition_classification;
condition_extract.artifact_times = [time_before_change,time_after_change];

cd(exp_path)

save('condition_extract.mat','condition_extract')

end

