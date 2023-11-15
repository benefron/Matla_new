function [allPer] = perExp(fileName,NoWhisking,numPer)
% Splits the tables based on the main experimntal conditions
load(fileName,"perTable")
mainTable = perTable(:,1:width(perTable)-2);
mainTable.Whisking = envelope(mainTable.Whisking);
alumTable = mainTable(strcmp("Aluminum",mainTable.("main labels")),:);
attenuatedTable = mainTable(strcmp("Attenuated",mainTable.("main labels")),:);
NoObjTable = mainTable(strcmp("No object",mainTable.("main labels")),:);
NoWhiskTable = mainTable(strcmp("No whisking",mainTable.("main labels")),:);

allPer.AlAt = updatedPer([attenuatedTable;alumTable],numPer,"main labels");
allPer.AlNO = updatedPer([NoObjTable;alumTable],numPer,"main labels");
allPer.AtNO = updatedPer([NoObjTable;attenuatedTable],numPer,"main labels");
% all experiments apart from FVB110_600 should run this part it is a
% control for the no whisking condition. and it doesn't compare the
% whisking in those conditions
if NoWhisking
    allPer.AlNW = updatedPer([NoWhiskTable(:,[1:end-2,end]);alumTable(:,[1:end-2,end])],numPer,"main labels");
    allPer.AtNW = updatedPer([NoWhiskTable(:,[1:end-2,end]);attenuatedTable(:,[1:end-2,end])],numPer,"main labels");
    allPer.NONW = updatedPer([NoWhiskTable(:,[1:end-2,end]);NoObjTable(:,[1:end-2,end])],numPer,"main labels");
    subTable = perTable(perTable.whiskingBinary == 0,[1:width(perTable)-4,width(perTable)-1]);
    subAlum= subTable(strcmp("NW_Aluminum",subTable.sublabels),:);
    subAttenuated= subTable(strcmp("NW_Attenuated",subTable.sublabels),:);
    subNoObject= subTable(strcmp("NW_No object",subTable.sublabels),:);
    allPer.noWhisking.AlAt = updatedPer([subAttenuated;subAlum],numPer,"sublabels");
    allPer.noWhisking.AlNO = updatedPer([subNoObject;subAlum],numPer,"sublabels");
    allPer.noWhisking.AtNO = updatedPer([subNoObject;subAttenuated],numPer,"sublabels");
end
saveFile = [fileName(1:10),'_permResults.mat'];
cd('Results/')
save(saveFile,"allPer")

end