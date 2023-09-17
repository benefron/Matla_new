function [binnedWhisking,whiskMatrics] = getWhiskingEnergy(experiment)
%This function accepts an object of class auditoryExp and extracts 1 sec
%bins of of whisking energy to perform statistical analysis across the
%conditions
    conditions = fieldnames(experiment.conditionVector);
    binnedWhisking.(conditions{1}) = getConditionWhisk(conditions{1});
    binnedWhisking.(conditions{2}) = getConditionWhisk(conditions{2});
    binnedWhisking.(conditions{3}) = getConditionWhisk(conditions{3});


    whiskMatrics.mean.(conditions{1}) = mean(binnedWhisking.(conditions{1}),2);
    whiskMatrics.mean.(conditions{2}) = mean(binnedWhisking.(conditions{2}),2);
    whiskMatrics.mean.(conditions{3}) = mean(binnedWhisking.(conditions{3}),2);

    whiskMatrics.std.(conditions{1}) = getSTDs(binnedWhisking.(conditions{1}));
    whiskMatrics.std.(conditions{2}) = getSTDs(binnedWhisking.(conditions{2}));
    whiskMatrics.std.(conditions{3}) = getSTDs(binnedWhisking.(conditions{3}));

    whiskMatrics.std.all = [whiskMatrics.std.(conditions{1});whiskMatrics.std.(conditions{2});whiskMatrics.std.(conditions{3})];
    whiskMatrics.mean.all = [whiskMatrics.mean.(conditions{1});whiskMatrics.mean.(conditions{2});whiskMatrics.mean.(conditions{3})];

    whiskMatrics.groups = [ones(1,length(whiskMatrics.std.(conditions{1}))),ones(1,length(whiskMatrics.std.(conditions{2})))*2,ones(1,length(whiskMatrics.std.(conditions{3})))*3];

    
    function binSTD = getSTDs(condMat)
        for i = 1:size(condMat,1)
            binSTD(i) = std(condMat(i,:));
        end
        binSTD = binSTD';
    end
    function [binned_whisking] = getConditionWhisk(condition)
        % Find the start and end of the whisking condition
        allDiffs = diff(experiment.conditionVector.(condition));
        values = [1,(find(allDiffs > 1)) + 1];
        values(2,:) = [find(allDiffs > 1),length(experiment.conditionVector.(condition))];

        % Find the corresponding camera frame from the whisking cam
        conditionEpochs = experiment.conditionVector.(condition)(values');
        conditionWhisking = conditionEpochs;
        kAll = []; 
        for ep = 1:length(conditionWhisking)
            [~,conditionWhisking(ep,1)] = min(abs(experiment.Cams.whisking.csv_aligned_frames-double(conditionEpochs(ep,1))));
            [~,conditionWhisking(ep,2)] = min(abs(experiment.Cams.whisking.csv_aligned_frames-double(conditionEpochs(ep,2))));
            kTemp = [conditionWhisking(ep,1):401:conditionWhisking(ep,2)];
            kTemp(end) = [];
            kAll = [kAll,kTemp];
        end

        kAll = kAll';
        binned_whisking = experiment.Whisking.Raw(repmat(0:399,length(kAll),1)+double(kAll));
    end

end