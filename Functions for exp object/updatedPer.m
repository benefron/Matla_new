function [pVal] = updatedPer(condTable,numPer,labelCoulmn)
meanComp = grpstats(condTable,labelCoulmn,'mean');
meanDiff = diff(meanComp(:,3:end));

permutedData = false(numPer, size(meanDiff, 2));

parfor i =1:numPer
    permutedLables = condTable;
    permutedLables.(labelCoulmn) = permutedLables.(labelCoulmn)(randperm(size(condTable,1)));
    meanPer = grpstats(permutedLables,labelCoulmn,'mean');
    meandiffPer = diff(meanPer(:,3:end));
    permutedData(i,:) = abs(table2array(meandiffPer)) >= abs(table2array(meanDiff));
end


pVal(1,:) = sum(permutedData)./numPer;

pVal(2,:) = sign(table2array(meanDiff));
end