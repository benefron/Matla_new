function [] = runModelOpt(fileName,modelType)
fullModel = load(fileName,modelType);
modelAll = fullModel.(modelType);

modelMain = modelAll(:,1:end-2);
subModel = modelAll(modelAll.whiskingBinary == 0,[1:end-3,end-1]);

cvMain = cvpartition(modelMain.("main labels"), 'Holdout', 0.2, 'Stratify', true);
MainModel.trainMain = modelMain(training(cvMain), :);
MainModel.testMain = modelMain(test(cvMain), :);

cvSub = cvpartition(subModel.sublabels, 'Holdout', 0.2, 'Stratify', true);
SubModel.trainSub = subModel(training(cvSub), :);
SubModel.testSub = subModel(test(cvSub), :);

t = templateTree('Reproducible',true); 

opts = struct('AcquisitionFunctionName', 'expected-improvement-plus', ...
              'MaxObjectiveEvaluations', 5, ...
              'UseParallel', true, ...
              'KFold', 10);
responseName = modelMain.Properties.VariableNames{end};
predictors = modelMain.Properties.VariableNames(1:end-1);
MainModel.mdl = fitcensemble(MainModel.trainMain(:, predictors), MainModel.trainMain.(responseName), ...
    'Method', 'LSBoost', ...
    'Learners', t, ...
    'OptimizeHyperparameters', 'auto', ...
    'HyperparameterOptimizationOptions', opts);
MainModel.yTest = MainModel.testMain.(responseName);
MainModel.yPred = MainModel.mdl.predict(MainModel.testMain(:,1:end-1));

responseName = subModel.Properties.VariableNames{end};
predictors = subModel.Properties.VariableNames(1:end-1);
SubModel.mdl = fitcensemble(SubModel.trainSub(:, predictors), SubModel.trainSub.(responseName), ...
    'Method', 'LSBoost', ...
    'Learners', t, ...
    'OptimizeHyperparameters', 'auto', ...
    'HyperparameterOptimizationOptions', opts);

SubModel.yTest = SubModel.testSub.(responseName);
SubModel.yPred = SubModel.mdl.predict(SubModel.testSub(:,1:end-1));

cd('Models/')
file_name = strjoin([fileName(1:10),modelType,'.mat'],'_')
save(file_name,"SubModel","MainModel")

end