function [unitTimesCon] = unitTimesCondition(experiment)
conVec = getConVector(experiment);
unitTimesCon.good.alum = getPerCon(experiment.Units.good.times,conVec.aluminum);
unitTimesCon.mua.alum = getPerCon(experiment.Units.mua.times,conVec.aluminum);

unitTimesCon.good.muted = getPerCon(experiment.Units.good.times,conVec.muted);
unitTimesCon.mua.muted = getPerCon(experiment.Units.mua.times,conVec.muted);

unitTimesCon.good.non = getPerCon(experiment.Units.good.times,conVec.non);
unitTimesCon.mua.non = getPerCon(experiment.Units.mua.times,conVec.non);

    function condMat = getPerCon(units,conVector)
        condMat.times = cell(size(units,1),1);
        isCon = ismember(units,conVector);
        condMat.FR = sum(isCon,2)/(length(conVector)/30000);
        for unit = 1:size(units)
            condMat.times{unit} = units(unit,isCon(unit,:));
        end
    end
end

