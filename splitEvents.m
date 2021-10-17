function [eventsCondition] = splitEvents(eventVector,conditionTimes)
% this function recieves a vector of times and matrix of all condition
% times and returns a vector of all times that are in range of the
% condition

eventsCondition = [];
%FVB103_757.experimentpredictedSound.events.aluminum = [];
for eve = 1:length(conditionTimes)
    tempEvents = eventVector(find(eventVector > conditionTimes(eve,1) & eventVector < conditionTimes(eve,2)));
    eventsCondition = [eventsCondition,tempEvents];
end
    
end

