% Unique levels of Mouse_ID and Cue_Type
mouse_ids = unique(data1.Mouse_ID);
cue_types = unique(data1.Cue_Type);
condition_types = unique(data1.Condition_Type);

% Define a colormap (one color for each mouse)
colormap = lines(length(mouse_ids));
colors = {colorBen.non,colorBen.aluminum,colorBen.muted}
% Create a new figure
figure;

mouseMarker = {'diamond','square','^'}
% For each Mouse_ID
for i = 1:length(mouse_ids)
    % For each Cue_Type
    for j = [1 2 3]
        % Filter the data
        df_mouse_cue = data1(data1.Mouse_ID == mouse_ids(i) & data1.Cue_Type == cue_types(j), :);
        
        % Calculate the mean outcome for each 'Condition_Type'
        interaction_means = groupsummary(df_mouse_cue, 'Condition_Type', 'mean', 'Outcome');
        reordredInteraction = [interaction_means.mean_Outcome(2),interaction_means.mean_Outcome(1),interaction_means.mean_Outcome(3)];
        % Plot
        plot(1:length(condition_types), reordredInteraction, 'Color', colors{j}+3*i*(0.05)*(1-colors{j}), 'DisplayName', ['Mouse ' char(mouse_ids(i)) ', Cue ' char(cue_types(j))],'Marker','*','LineStyle',':','MarkerFaceColor',colors{j}+3*i*(0.05)*(1-colors{j}),'MarkerSize',5);
        hold on;
    end
end
% Calculate and plot the overall mean outcome for each 'Condition_Type' and 'Cue_Type'
for j = [3 2 1]
    % Filter the data
    df_cue = data1(data1.Cue_Type == cue_types(j), :);
    
    % Calculate the mean outcome for each 'Condition_Type'
    overall_means = groupsummary(df_cue, 'Condition_Type', 'mean', 'Outcome');
    overall_std = groupsummary(df_cue, 'Condition_Type', 'std', 'Outcome');
    reordredMeans= [overall_means.mean_Outcome(2),overall_means.mean_Outcome(1),overall_means.mean_Outcome(3)];
    reordredSTD= [overall_std.std_Outcome(2),overall_std.std_Outcome(1),overall_std.std_Outcome(3)];
    overall_se = reordredSTD/sqrt(3);
    % Plot
    %boxplot(df_cue.Outcome,df_cue.Condition_Type)
    errorbar(reordredMeans,overall_se,'DisplayName', ['Overall, Cue ' char(cue_types(j))],'LineStyle','none','Marker','o','MarkerFaceColor',colors{j},'Color',colors{j},'MarkerSize',6)
    %plot(1:length(condition_types), overall_means.mean_Outcome,'DisplayName', ['Overall, Cue ' char(cue_types(j))],'LineStyle','none','Marker','.','Color',colors{j},MarkerSize=15);
    hold on;
end

% Legend
legend({'Individual mouse','','','','','','','','','Attenuated (no-go)','Aluminum(go)','No object(no go)'});

% Labels
% xlabel('Condition Type');
xticks(1:length(condition_types));
xticklabels({'Task','Catch trials','White noise'});
ylabel('Lick rate');
% title("Interaction between Cue Type and Condition Type");
xlim([0.8 3.2])
ylim([0 1.1])
% Show the plot
hold off;
set(gca,'FontSize',12)
set(gcf,'Position',[2267,71,1058,807]);
set(gca,'FontName','SansSerif')
