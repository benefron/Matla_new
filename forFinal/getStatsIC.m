function [pVal] = getStatsIC(rasterStruct,afterSec)
% get's a structre of the rasters of all units and preforms a paired
% statistical test between base line and time after galvo movement
windowBefore = [1:length(rasterStruct{1,1})/2-1];
beforeSec = length(windowBefore)/30000;
afterWindow = [(length(rasterStruct{1,1})/2)+1:(length(rasterStruct{1,1})/2+1)+afterSec*30000-1];
pVal = [];


for unit = 1:size(rasterStruct,1)
    alumTemp = rasterStruct{unit,1};
    mutTemp = rasterStruct{unit,2};
    nonTemp = rasterStruct{unit,3};
    % baselinealum = sum(alumTemp(:,windowBefore),2)/beforeSec;
    % baselinemut = sum(mutTemp(:,windowBefore),2)/beforeSec;
    % baselinenon = sum(nonTemp(:,windowBefore),2)/beforeSec;
    % meanBase = mean(baselineFR);
    alumVec = sum(alumTemp(:,afterWindow),2)/afterSec;
    mutVec = sum(mutTemp(:,afterWindow),2)/afterSec;
    nonVec = sum(nonTemp(:,afterWindow),2)/afterSec;
    forStat = [alumVec;mutVec;nonVec];
    groups = [repmat({'Aluminum'},length(alumVec),1);repmat({'Attenuated'},length(mutVec),1);repmat({'No Object'},length(nonVec),1)];
    [p,tbl,stats] = kruskalwallis(forStat,groups,'off');
    [mulCo,h] = multcompare(stats,'Display','off');
    pVal(unit,1) = p;
    pVal(unit,2:4) = mulCo(:,6)';
end

end