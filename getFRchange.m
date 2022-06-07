function [FR,stat] = getFRchange(experiment,binSize)
goodUnits = [1:length(experiment.Units.perCondition.good.non.times)];
muaUnits = [1:length(experiment.Units.perCondition.mua.non.times)]+goodUnits(end);
all_alum = [experiment.Units.perCondition.good.alum.times;experiment.Units.perCondition.mua.alum.times];
all_mut = [experiment.Units.perCondition.good.muted.times;experiment.Units.perCondition.mua.muted.times];
all_non = [experiment.Units.perCondition.good.non.times;experiment.Units.perCondition.mua.non.times];

FR.Aluminum = getCondFR(experiment.conditionVector.aluminum,all_alum,binSize);
FR.Muted = getCondFR(experiment.conditionVector.muted,all_mut,binSize);
FR.Non = getCondFR(experiment.conditionVector.non,all_non,binSize);

for test = 1:length(all_alum)
    [stat.Right(test,1),stat.Right(test,2)] = ranksum(FR.Aluminum.Values(test,:),FR.Non.Values(test,:),'tail','right','alpha',0.025);
    [stat.Left(test,1),stat.Left(test,2)] = ranksum(FR.Aluminum.Values(test,:),FR.Non.Values(test,:),'tail','left','alpha',0.025);
end


    function condFR = getCondFR(condition,units,binSize)
        binSize = binSize*30000;
        timeIndex(:,2) = find(diff(condition) > 2 * binSize);
        timeIndex(:,1) = [1;timeIndex(1:end-1,2)+1];
        binBorders = [];
        time2discard = [];
        for epoch = 1:length(timeIndex)
           tempBin = [timeIndex(epoch,1):binSize:timeIndex(epoch,2)];
           time2discard = [time2discard,tempBin(end):timeIndex(epoch,2)];
           binBorders = [binBorders,tempBin];
        end
        time2discard = condition(time2discard);
        binBorders = condition(binBorders);
        condFR.Values = zeros(length(units),length(binBorders)-1);
        for un = 1:length(units)
            unit = units{un};
            unitsOut = ~ismember(unit,time2discard);
            unit = unit(unitsOut);
            h = histogram(unit,'BinEdges',binBorders);
            condFR.Values(un,:) = h.Values;
        end
        condFR.Bins = binBorders;
    end

end

