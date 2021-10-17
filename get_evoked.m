function [evoked] = get_evoked(varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
funCase = 3;
run(varargin{1});
switch openEphDataType
case 1 % openephys format
        dataSPK = load_open_ephys_data_faster(speakerChannel)';
case 2 % .dat bin file format
        dataSPK = varargin{2}.mapped(speakerChannel,1:end);
    case 3
        dataSPK = speakerChannel';
end

dataSPK = dataSPK - mean(dataSPK);

all_aboveSPK = find(abs(dataSPK) > ThrSPK); % Finds all the frames were the arrduino sent a 5v signal 
all_belowSPK = find(abs(dataSPK) < ThrSPK); % Finds all events were arduino did not send a 5V signal
all_startsSPK = intersect(all_aboveSPK',(all_belowSPK+1)'); % Finds the first frame of every 5V signal
all_startsSPK(all_startsSPK > evokedTime(2) * SF | all_startsSPK < evokedTime(1) * SF) = [];


wind_mat = (ones(length(all_startsSPK),0.1*SF) .* [-0.1*SF:1:-1]) + all_startsSPK; % create a clean matrix to find the values before the whisk_onsey
    
wind_dat_onset = abs(dataSPK(wind_mat)); % find the value for all frames before the onset of whisking in absolute value 
[rows,~] = find(wind_dat_onset > ThrSPK); rows = unique(rows);% finds the rows where the frames before had a values above threshold
all_startsSPK(rows) = [];
evoked = all_startsSPK;


to_plot = 0;
switch to_plot
    case 0
    case 1
        figure; plot (dataSPK);
        hold on
        for i=1:length(evoked)
            xline(evoked(i),'r');
        end
end



end

