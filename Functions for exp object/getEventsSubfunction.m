function [eventsCond,eventsEp,totalTime,eventsWaveform] = getEventsSubfunction(dataPascal,epochs,th,aviTimes)
eventsCond.real = [];
eventsEp.real = [];
eventsCond.shuffled = [];
eventsEp.shuffled = [];
totalTime = 0;
eventsWaveform = {};
k = 1;
for ep  = 1:sum(~isnan(epochs(:,2)))
    endIndex = min(length(dataPascal), epochs(ep,2) + 500);
    smoothEpoch = dataPascal(epochs(ep,1)-500:endIndex);
    endTrim = 500 - max(0, (epochs(ep,2) + 500) - length(dataPascal));
    smoothEpoch = envelope(smoothEpoch,325,'peak');
    smoothEpoch = smoothEpoch(501:end-endTrim);
    %eventsEP = abs(smoothEpoch) > th;
    events_diff = diff([0; eventsEP; 0]);
    startsPre = find(events_diff == 1);
    endsPre = find(events_diff == -1) - 1;
    shortEvents = endsPre-startsPre < 500;
    startsPre(shortEvents) = [];
    endsPre(shortEvents) = [];
    [starts,ends] = mergeCloseEvents(startsPre, endsPre, 250);
    % if  ep == 1
    %     figure;
    %     plot(dataPascal(epochs(ep,1):epochs(ep,2)));
    %     try
    %         for i = 1:30
    %             xline(starts(i),'r')
    %             xline(ends(i),'g')
    %             hold on;
    %         end
    %     catch
    %         for i = 1:length(starts)
    %             xline(starts(i),'r')
    %             xline(ends(i),'g')
    %             hold on;
    %         end
    %     end
    %     plot(smoothEpoch)
    % end
    conditionV = epochs(ep,1):epochs(ep,2);
    startsRand = randperm(length(conditionV),length(starts));
    starts = conditionV(starts);
    ends = conditionV(ends);
    startsEp = round(aviTimes(starts));
    startsRand = conditionV(startsRand);
    startsEpRand = round(aviTimes(startsRand));
    eventsCond.real = [eventsCond.real,starts];
    eventsEp.real = [eventsEp.real,startsEp];
    eventsCond.shuffled = [eventsCond.shuffled,startsRand];
    eventsEp.shuffled = [eventsEp.shuffled,startsEpRand];
    totalTime = totalTime + length(conditionV)/250000;
    for i = 1:length(starts)
        eventsWaveform{k} = [starts(i):ends(i)];
        k=k+1;
    end

end

    function [mergedStarts, mergedEnds] = mergeCloseEvents(starts, ends, minGap)
        % Initialize merged events as the original events
        mergedStarts = starts;
        mergedEnds = ends;

        % Calculate gaps between consecutive events
        gaps = [Inf; starts(2:end) - ends(1:end-1)]; % Prepend Inf to keep the first event

        % Find events to merge based on gap criteria
        mergeIndices = find(gaps <= minGap & gaps > 0);

        % Adjust ends to extend to the next event's end, removing merged starts and ends
        while ~isempty(mergeIndices)
            % Adjust the end of the current event to the end of the next one
            mergedEnds(mergeIndices - 1) = mergedEnds(mergeIndices);

            % Remove the merged events' starts and subsequent redundant ends
            mergedStarts(mergeIndices) = [];
            mergedEnds(mergeIndices - 1) = []; % Keep the updated end, remove the next

            % Recalculate gaps to find new events to merge
            gaps = [Inf; mergedStarts(2:end) - mergedEnds(1:end-1)];
            mergeIndices = find(gaps <= minGap & gaps > 0);
        end
    end
end


