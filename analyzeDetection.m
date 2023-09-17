function [results] = analyzeDetection(events_vec,id)
% This function takes the vector of events from the test session of the
% whisking sounds detection task and divides it into the different controls
% and catagories to calculate hit rate, false alarm rate and d'
just_relevant = events_vec(logical(sum((events_vec == [444,404,600,601,0,1,2,3,4,5])')));
TrialType = [404,444];
only400 = ismember(just_relevant,TrialType);
ind400 = find(only400==1);
justReshape = [];
% diff400 = diff(ind400);
% toDel = find(diff400==1);
% toDel = ind400(toDel);
% just_relevant(toDel) = [];
k=1;
for t = 1:length(ind400)
    try
        if and(or(just_relevant(ind400(t)+1) == 600,just_relevant(ind400(t)+1) == 601),just_relevant(ind400(t)+2) < 6)
            justReshape(1,k) = just_relevant(ind400(t));
            justReshape(2,k) = just_relevant(ind400(t)+1);
            justReshape(3,k) = just_relevant(ind400(t)+2);
            k = k+1;
        end
    catch
    end
end

%justReshape = reshape(just_relevant,3,length(just_relevant)/3);
white_noise = justReshape(:,find(justReshape(1,:)==444));
nonWhite  = justReshape(:,find(justReshape(1,:)==404));
catchTrials = nonWhite(3,find(nonWhite(2,:)==601));
regularTrials = nonWhite(3,find(nonWhite(2,:)==600));
results.whiteNoise = getRates(white_noise(3,:));
results.regularTrials = getRates(regularTrials);
results.catchTrials = getRates(catchTrials);
% create empty vectors for table
mouse_id = repmat({id},length(justReshape),1);
cue_type = cell(length(justReshape),1);
Outcome = zeros(length(justReshape),1);
conditionType = cell(length(justReshape),1);

% allocate the conditions and cue to the table
missT = justReshape(3,:) == 0;
cue_type(missT) = {'Aluminum Go'};
Outcome(missT) = 0;
hitT = justReshape(3,:) == 1;
cue_type(hitT) = {'Aluminum Go'};
Outcome(hitT) = 1;
CRT = justReshape(3,:) == 2;
cue_type(CRT) = {'No Go'};
Outcome(CRT) = 0;
FAT = justReshape(3,:) == 3;
cue_type(FAT) = {'No Go'};
Outcome(FAT) = 1;
AtCRT = justReshape(3,:) == 4;
cue_type(AtCRT) = {'Attenuated No Go'};
Outcome(AtCRT) = 0;
AtFAT = justReshape(3,:) == 5;
cue_type(AtFAT) = {'Attenuated No Go'};
Outcome(AtFAT) = 1;

whiteNoiseInd = find(justReshape(1,:)==444);
catchTrailInd = find(justReshape(2,:)==601);

conditionType(whiteNoiseInd) = {"White noise"};
conditionType(catchTrailInd) = {"Catch Trial"};
conditionType(cellfun(@isempty,conditionType)) = {"Regular"};



exp_table = table(mouse_id, cue_type, conditionType, Outcome, 'VariableNames', {'Mouse_ID', 'Cue_Type', 'Condition_Type', 'Outcome'});
results.table = exp_table;


    function [condStats] = getRates(outcomeVec)
        go.hits = sum(outcomeVec == 1);
        go.miss = sum(outcomeVec == 0);
        noGo.CR = sum(outcomeVec == 2);
        noGo.FA = sum(outcomeVec == 3);
        atNoGo.CR = sum(outcomeVec == 4);
        atNoGo.FA = sum(outcomeVec == 5);
        go.hitRate = go.hits/(go.hits+go.miss);
        noGo.FARate = noGo.FA/(noGo.FA + noGo.CR);
        atNoGo.FARate = atNoGo.FA/(atNoGo.FA + atNoGo.CR);
        noGo.CRRate = noGo.CR/(noGo.FA + noGo.CR);
        atNoGo.CRRate = atNoGo.CR/(atNoGo.FA + atNoGo.CR);
        dPrime.goNoGo = norminv(go.hitRate) - norminv(noGo.FARate);
        dPrime.goAttNoGo = norminv(go.hitRate) - norminv(atNoGo.FARate);
        condStats.go = go;
        condStats.atNoGo = atNoGo;
        condStats.noGo = noGo;
        condStats.dPrime = dPrime;
    end

end