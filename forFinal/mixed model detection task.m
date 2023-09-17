% Load the data
data = readtable('/media/ben/Extreme SSD/Detection task/F11134.csv');

% Convert categorical variables to categories
data.Cue_Type = categorical(data.Cue_Type);
data.Condition_Type = categorical(data.Condition_Type);
data.Mouse_ID = categorical(data.Mouse_ID);
data1 = data;
% Set 'No Go' as the reference level for 'Cue_Type'
data1.Cue_Type = reordercats(data.Cue_Type, {'No Go','Aluminum Go','Attenuated No Go'});
data1.Condition_Type = reordercats(data.Condition_Type, {'Catch Trial','Regular','White noise'});

% Define the model
formula = 'Outcome ~ Cue_Type * Condition_Type + (1|Mouse_ID)';
glme = fitglme(data1, formula, 'Distribution', 'Binomial', 'Link', 'logit', 'FitMethod', 'Laplace');

% Get the fixed effects
[~,~,FEParams] = fixedEffects(glme);
disp(glme)
exp(FEParams.Estimate)
% Calculate probabilities
probabilities = 1 ./ (1 + exp(-FEParams.Estimate));

% Get the standard errors
stderr = FEParams.SE;
labels_all = {'Intercept(No object(No go)/Catch)','Aluminum(Go)/catch','Attenuated(No go)/catch','White noise/No object(No go)','Task/No object(No go)','Aluminum X White noise','Attenuated X White noise','Aluminum X Task','Attenuated X Task'}
% Plot the probabilities with error bars representing the standard errors
figure;

errorbar(1:numel(FEParams.Estimate), FEParams.Estimate, FEParams.Lower,FEParams.Upper,'.','MarkerSize',15);
%xlabel('Variable');
ylabel('log odd');
title('ref: No object (No go),catch');
xticks(1:numel(FEParams.Estimate));
xticklabels(labels_all);
xtickangle(90);
xlim([0 10])
yline(0,'--','Color','r')

% Display the model summary
disp(glme)

