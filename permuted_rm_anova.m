function [p_value, posthoc_results] = permuted_rm_anova(filename, n_permutations)
    % Load data
    data = readtable(filename);
    
    % Convert data to wide format
    data_wide = unstack(data, 'Outcome', 'Condition_Type', 'GroupingVariables', {'Mouse_ID', 'Cue_Type'});
    
    % Construct formula for fitrm based on the number of Condition_Type levels
    conditionLevels = unique(data.Condition_Type);
    formula = sprintf('Outcome_%s-Outcome_%s~Cue_Type', conditionLevels{:});
    
    % Fit a repeated measures model
    rm = fitrm(data_wide, formula, 'WithinDesign', table(conditionLevels));
    
    % ... (rest of the function remains unchanged)

    % Compute observed F-statistic for main test
    ra = ranova(rm);
    observed_F = ra.F(1);
    
    % Permutation test for main effect
    F_values = zeros(n_permutations, 1);
    for i = 1:n_permutations
        shuffled_data = data;
        for mouse = unique(data.Mouse_ID)'
            idx = data.Mouse_ID == mouse;
            shuffled_data.Outcome(idx) = data.Outcome(randperm(sum(idx)));
        end
        rm_permuted = fitrm(shuffled_data, 'Outcome~Cue_Type*Condition_Type', 'WithinDesign', table(shuffled_data.Cue_Type, shuffled_data.Condition_Type), 'WithinModel', 'Cue_Type*Condition_Type');
        ra_permuted = ranova(rm_permuted);
        F_values(i) = ra_permuted.F(1);
    end
    
    % Compute p-value for main effect
    p_value = sum(F_values >= observed_F) / n_permutations;
    
    % Permutation-based post-hoc pairwise comparisons
    groups = unique(data.Cue_Type);  % or data.Condition_Type, depending on which factor you're testing
    n_groups = length(groups);
    posthoc_results = table;
    for i = 1:n_groups
        for j = (i+1):n_groups
            group1 = data(data.Cue_Type == groups(i), :);  % or data.Condition_Type
            group2 = data(data.Cue_Type == groups(j), :);  % or data.Condition_Type
            observed_diff = mean(group1.Outcome) - mean(group2.Outcome);
            
            diff_values = zeros(n_permutations, 1);
            for k = 1:n_permutations
                permuted_group1 = group1;
                permuted_group2 = group2;
                for mouse = unique(data.Mouse_ID)'
                    idx1 = group1.Mouse_ID == mouse;
                    idx2 = group2.Mouse_ID == mouse;
                    combined = [group1.Outcome(idx1); group2.Outcome(idx2)];
                    shuffled = combined(randperm(length(combined)));
                    permuted_group1.Outcome(idx1) = shuffled(1:sum(idx1));
                    permuted_group2.Outcome(idx2) = shuffled((sum(idx1)+1):end);
                end
                diff_values(k) = mean(permuted_group1.Outcome) - mean(permuted_group2.Outcome);
            end
            
            p = sum(diff_values >= observed_diff) / n_permutations;
            posthoc_results = [posthoc_results; table(groups(i), groups(j), observed_diff, p)];
        end
    end
    
    % Bonferroni correction
    posthoc_results.p = posthoc_results.p * size(posthoc_results, 1);
    posthoc_results.p(posthoc_results.p > 1) = 1;
end
