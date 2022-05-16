function conditionVector = getConVector(experiment)
artifactStart = experiment.Conditions.artifact_times(2) * 30;
artifactEnd = experiment.Conditions.artifact_times(1) * 30;
all_times = [experiment.Conditions.all_changes + artifactStart,[experiment.Conditions.all_changes(2:end);length(experiment.sound_events.fullSound)]-artifactEnd];
aluminumTimes = all_times(experiment.Conditions.condition_classification == 1,:);
mutedTimes = all_times(experiment.Conditions.condition_classification == 2,:);
nonTimes = all_times(experiment.Conditions.condition_classification == 3,:);
% create idx vector for aluminum
conditionVector.aluminum = [];
for al = 1:length(aluminumTimes)
    tempInd = [aluminumTimes(al,1):aluminumTimes(al,2)];
    conditionVector.aluminum = [conditionVector.aluminum,tempInd];
end
conditionVector.muted = [];
for al = 1:length(mutedTimes)
    tempInd = [mutedTimes(al,1):mutedTimes(al,2)];
    conditionVector.muted = [conditionVector.muted,tempInd];
end
conditionVector.non = [];
for al = 1:length(nonTimes)
    tempInd = [nonTimes(al,1):nonTimes(al,2)];
    conditionVector.non = [conditionVector.non,tempInd];
end
end

