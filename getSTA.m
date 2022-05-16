function [STA] = getSTA(experiment,timeWindow,avgVector,vecType)
all_alum = [experiment.Units.perCondition.good.alum.times;experiment.Units.perCondition.mua.alum.times];
all_mut = [experiment.Units.perCondition.good.muted.times;experiment.Units.perCondition.mua.muted.times];
all_non = [experiment.Units.perCondition.good.non.times;experiment.Units.perCondition.mua.non.times];
goodUnits = [1:length(experiment.Units.perCondition.good.non.times)];
muaUnits = [1:length(experiment.Units.perCondition.mua.non.times)]+goodUnits(end);
if vecType == 1
    sf = 400;
else
    sf = 30000;
end

STA.Aluminum = getSTACond(all_alum);
STA.Muted = getSTACond(all_mut);
STA.Non = getSTACond(all_non);
STA.goodU = goodUnits;
STA.muaU = muaUnits;


    function [unitSTA] = getUnitSTA(spkTrain,timeWindow)
        
        timeWhisk = single([1:timeWindow*sf] - (timeWindow*sf)/2);
        unitWhisk = single(ones(length(spkTrain),length(timeWhisk)).*timeWhisk);
        if vecType == 1
            spkTrainWhisk = single(zeros(length(spkTrain),1));
            for i = 1:length(spkTrain)
                [~,spkTrainWhisk(i)] = min(abs(experiment.Cams.whsiking.csv_aligned_frames-spkTrain(i)));
            end
            unitWhisk = (timeWhisk+spkTrainWhisk);
        else
            unitWhisk = (timeWhisk+spkTrain');
            
        end
        unitSTA = mean(avgVector(unitWhisk));
    end
    function [condSTA] = getSTACond(conditionSPK)
        condSTA = zeros(length(conditionSPK),timeWindow*sf);
        for unit = 1:length(conditionSPK)
           tempCond = getUnitSTA(conditionSPK{unit},timeWindow); 
           condSTA(unit,:) = tempCond;
        end
        
    end



end
