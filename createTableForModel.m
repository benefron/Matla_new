function [ClassificationTable] = createTableForModel(expID)
%Function for aliginig sound events and camera frames and extracting main
%features for model fitting
%   Input recives exp ID (animalID_recordingDepth) to upload the following: 
%           - experiment object of the class 'auditory_ephys_exp'
%           - Sound events labels from .txt file
%           - Analysis of movement energy from faceMap app .mat
%   The function output a formatted table with the following fields (27) for each event:
%           - Total energy vector in a 60ms time window (for 2 ROIs)
%           - Mean total motion energy (for 2 ROIs)
%           - STD of total motion energy (for 2 ROIs)
%           - 10 first motion masks vector in a 60ms time window (for 2
%             ROIs)
%           - Classification


%% load all relevant information
cd('/home/ben/Z')
load(['/home/ben/Z/Experiment objects/',expID, '.mat']);
labeled_events = openAVisoftLabels(['/home/ben/Z/Labels/',expID,'.txt']);
load(['/home/ben/Z/All_whiskers_videos/Analyzed/',expID,'.mat'],'motion_1','motion_2','motSVD_1','motSVD_2');
AvisoftTimes = temp_exp.Avisoft_times;
CamTimes = temp_exp.Cams.Whisking;
clear temp_exp

%% get silence times equal to alum times
alum_tol = uniquetol(labeled_events.Alu,0.2/max(labeled_events.Alu));
eventQ = length(alum_tol);
all_silence = zeros(1,(length(labeled_events.Silence) * 500));
k = [1,500];
for epoch = 1:length(labeled_events.Silence)
    a = labeled_events.Silence(epoch,1);
    b = labeled_events.Silence(epoch,2);
    all_silence(k(1):k(2)) = (b-a).*rand(500,1) + a;
    k = k + 500;
end
labeled_events.SilenceRandom = all_silence(randi(length(all_silence),eventQ,1))';

%% get the labeled data in aligned samples
all_classification = [ones(eventQ,1);zeros(eventQ,1)];
all_events = [alum_tol;labeled_events.SilenceRandom];
all_events_sample = round(all_events * AvisoftTimes.SF);
all_events_aligned = AvisoftTimes.alignedFrames(all_events_sample);
all_events_openephys_sample = round(all_events_aligned/AvisoftTimes.sampling_ratio);

%% find the frame closest to the aligned samples
for event = 1:eventQ*2
    [~,closestIndex(event)] = min(abs(CamTimes.csv_aligned_frames - all_events_openephys_sample(event)));
end


%% create the features for model
timeVector = [-0.01 * CamTimes.fps : 0.05 * CamTimes.fps-1];
ClassificationTable = table;
motion_1_labeled = motion_1((ones(event,length(timeVector)).*timeVector) + closestIndex');
motion_2_labeled = motion_2((ones(event,length(timeVector)).*timeVector) + closestIndex');
ClassificationTable.motion_1_mean = mean(motion_1_labeled,2);
ClassificationTable.motion_2_mean = mean(motion_2_labeled,2);
ClassificationTable.motion_1_STD = std(motion_1_labeled,0,2);
ClassificationTable.motion_2_STD = std(motion_2_labeled,0,2);
if CamTimes.fps == 300
    motion_1_labeled = single(resample(double(motion_1_labeled),400,300,'Dimension',2));
    motion_2_labeled = single(resample(double(motion_2_labeled),400,300,'Dimension',2));
end
ClassificationTable.motion_1 = motion_1_labeled;
ClassificationTable.motion_2 = motion_2_labeled;

for mot = 1:10
    tempMot_1 = motSVD_1(:,mot);
    tempMot_2 = motSVD_2(:,mot);
    tempMot_1_resample = tempMot_1((ones(event,length(timeVector)).*timeVector) + closestIndex');
    tempMot_2_resample = tempMot_2((ones(event,length(timeVector)).*timeVector) + closestIndex');
    if CamTimes.fps == 300
        tempMot_1_resample = single(resample(double(tempMot_1_resample),400,300,'Dimension',2));
        tempMot_2_resample = single(resample(double(tempMot_2_resample),400,300,'Dimension',2));     
    end
    
    ClassificationTable.(['MotSVD_1_',num2str(mot)]) = tempMot_1_resample;
    ClassificationTable.(['MotSVD_2_',num2str(mot)]) = tempMot_2_resample;
end

ClassificationTable.Class = all_classification;
end

