function [predictionTable,zeroTime] = createTableForPrediction(expID)
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
load(['Experiment objects/',expID, '.mat']);
%labeled_events = openAVisoftLabels(['Labels/',expID,'.txt']);
load(['All_whiskers_videos/Analyzed/',expID,'.mat'],'motion_1','motion_2','motSVD_1','motSVD_2');
CamTimes = temp_exp.Cams.Whisking;
clear temp_exp

%% 
timeAfter = CamTimes.fps * 0.05-1;
totalFrameWindow = 1:CamTimes.fps * 0.06;
timeMatrix = (CamTimes.csv_start_frame:length(CamTimes.csv_aligned_frames) - timeAfter)';
timeMatrix = totalFrameWindow + timeMatrix-1;
%% create the features for model
predictionTable = table;
motion_1_labeled = motion_1(timeMatrix);
motion_2_labeled = motion_2(timeMatrix);
predictionTable.motion_1_mean = mean(motion_1_labeled,2);
predictionTable.motion_2_mean = mean(motion_2_labeled,2);
predictionTable.motion_1_STD = std(motion_1_labeled,0,2);
predictionTable.motion_2_STD = std(motion_2_labeled,0,2);
if CamTimes.fps == 300
    motion_1_labeled = single(resample(double(motion_1_labeled),400,300,'Dimension',2));
    motion_2_labeled = single(resample(double(motion_2_labeled),400,300,'Dimension',2));
end
predictionTable.motion_1 = motion_1_labeled;
predictionTable.motion_2 = motion_2_labeled;

for mot = 1:10
    tempMot_1 = motSVD_1(:,mot);
    tempMot_2 = motSVD_2(:,mot);
    tempMot_1_resample = tempMot_1(timeMatrix);
    tempMot_2_resample = tempMot_2(timeMatrix);
    if CamTimes.fps == 300
        tempMot_1_resample = single(resample(double(tempMot_1_resample),400,300,'Dimension',2));
        tempMot_2_resample = single(resample(double(tempMot_2_resample),400,300,'Dimension',2));     
    end
    
    predictionTable.(['MotSVD_1_',num2str(mot)]) = tempMot_1_resample;
    predictionTable.(['MotSVD_2_',num2str(mot)]) = tempMot_2_resample;
end
zeroTime = timeMatrix(:,5);
zeroTime(:,2) = CamTimes.csv_aligned_frames(zeroTime);
%predictionTable.Class = all_classification;
end

