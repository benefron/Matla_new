function [eventsMatrix] = meanMatrix(events,timeWindow,sf)
% This functions outputs a matrix with each row a vector of samples around
% an event
timeVector = [-0.5*timeWindow*sf+1:0.5*timeWindow*sf];
emptyMatrix = ones(length(events),length(timeVector));
eventsMatrix = emptyMatrix .*timeVector;
try
    eventsMatrix = eventsMatrix + events;
catch
    eventsMatrix = eventsMatrix + events';
end
end

