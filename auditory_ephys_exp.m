classdef auditory_ephys_exp
    % Processed data from ephys experiment for auditory-tectile project
    % starting Dec 2020
    %   This class contains the unit data together with the relevent analog
    %   channels and timestamps of digital events. 
    %   The class contains the following info:
    %       Pathways to data files
    %       Graph colors
    %       Unit times and ID
    %       Artifact times
    %       data from params_file
    %       Synced camera frame times
    %       Condition times
    %       Event times
    %       audio event classification and times
    %       Predicted audio times from model
    
    properties
        experiment_ID % name and identifyer for experiment
        Pathways    % all relevant pathways to the data
        digital_channels    % Digital channel number from open ephys
        analog_channels     % analog channel numbers form open ephys
        artifacts_parameters    % parameters for the time around artifact to clear data
        experiment_metadata     % data about the animal: ID, DOB,Recording depth, date of experiment
        Units   % Units times and identity !!! align with avi if needed
        Conditions  % Times of different experiment conditions
        Cams    % Timestamps of camrea frames !!!align with avi if needed
        SF % experiment sampling rate
        color % The RGB colors for the figures
        sound_events % event times extracted from sound and classification
        PSTHs
        sound
        running
        evokedEvents;
        conditionVector;
        Whisking
        timeStamp
        analysisType
        time2discard
        SDF
        
        
    end
    
    methods 
        function obj = auditory_ephys_exp(parameters_file) 
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            run(parameters_file);
            obj.experiment_ID = expName;
            obj.Pathways.exp = exp_path;
            obj.Pathways.parameters = params_path;
            obj.Pathways.sorting = sorting_path;
            obj.Pathways.auditory = auditory_file_path;
            obj.Pathways.main = folder_path;
            obj.Pathways.FPGA = FPGA_folder;
            obj.SF = SF;
            obj.color.Aluminum = [0.0039 * 210 , 0.0039 * 30 , 0.0039 * 75]; % D21E4B adobe color wheel
            obj.color.noObject = [0.0039 * 232 , 0.0039 * 134 , 0.0039 * 65]; % 284B59 adobe color wheel
            obj.color.Muted = [0.0039 * 40, 0.0039 * 75, 0.0039 * 89]; % E88641 adobe color wheel
            % all data channels from open ephys
            obj.digital_channels.arduino = arduino_channel;
            obj.digital_channels.face_cam = face_camera;
            obj.digital_channels.whisker_cam = Whiskers_cam;
            obj.analog_channels.encoder = encoder;
            obj.analog_channels.mic = mic;
            obj.analog_channels.speaker = speaker;
            obj.artifacts_parameters.time_before = time_before_change;
            obj.artifacts_parameters.time_after = time_after_change;
            obj.artifacts_parameters.lengthOnOff = length_on_off;
            obj.experiment_metadata.animalID = Animal_ID;
            obj.experiment_metadata.DOB = DOB;
            obj.experiment_metadata.expDate = Date_of_exp;
            obj.experiment_metadata.recording_num = recording_num;
            try
                obj.time2discard = time_to_discard;
            catch
            end
            %obj.Conditions = obj.getConditionsTimes;
            %obj.Units = sorted_from_phy(obj.Pathways.sorting);
            %obj.Cams.whsiking = obj.syncCamera(1);
            %obj.Cams.face = obj.syncCamera(2);
            %obj.evokedEvents = obj.getEvoked;
            %obj.Units.audUnits = find(obj.evokedEvents.pVal(:,2)==1);
        end
        
        function condition_extract = getConditionsTimes(obj)
            cd([obj.Pathways.exp,'/experiment',num2str(obj.experiment_metadata.recording_num),'/recording1/events/',obj.Pathways.FPGA,'TTL_1/'])
            open_ephys_TTL = load_open_ephys_binary([obj.Pathways.exp,'/experiment',num2str(obj.experiment_metadata.recording_num),'/recording1/structure.oebin'],'events',1);
            channel_states = open_ephys_TTL.Data;
            timestamps = open_ephys_TTL.Timestamps;
            cd([obj.Pathways.exp,'/experiment',num2str(obj.experiment_metadata.recording_num),'/recording1/continuous/',obj.Pathways.FPGA])
            all_timestamps = readNPY('timestamps.npy');
            
            timestamps = timestamps - all_timestamps(1);
            OnTimestamps = timestamps(channel_states == obj.digital_channels.arduino);
            
            isLastOn = [(OnTimestamps(2:end) - OnTimestamps(1:end-1)) > ((obj.artifacts_parameters.lengthOnOff * 2.5) * (obj.SF/1000));1];
            positions = find(isLastOn);
            condition_classification = [positions(1);(positions(2:end) - positions(1:end-1))];
            all_changes = OnTimestamps(positions);
            condition_classification(diff(all_changes) < 850000) = [];
            all_changes(diff(all_changes) < 850000) = [];
            k = 1;
%             for i = 1:length(all_changes)
%                startTime = all_changes(i);
%                if i==length(all_changes)
%                    endTime = length(all_timestamps);
%                else
%                    endTime = all_changes(i+1);
%                    
%                end
%                if endTime-startTime > 1200000
%                    f = warndlg(["there is a missing change of condition after condition: ",num2str(i)]);
%                    condition_extract.missed(k) = i;
%                    k = k+1;
%                    
%                end
%                    
%             end
            condition_extract.all_changes = all_changes;
            condition_extract.condition_classification = condition_classification;
            condition_extract.artifact_times = [obj.artifacts_parameters.time_before,obj.artifacts_parameters.time_after];
            
            cd(obj.Pathways.exp)

            save('condition_extract.mat','condition_extract')
        end
        
        function samples_aluminum = identifyEvents(obj)
            % finds the self generated events from the aluminum condition
            [data,fs] = audioread(obj.Pathways.auditory,'native');% load audio synced as native
            condition_times = double(obj.Conditions.all_changes)/obj.SF;
            condition_times = round(condition_times*fs);
            condition_times(:,2) = [condition_times(2:end,1);length(data)]; % Get's the begining and end of each condition changes
            alum_times = condition_times(obj.Conditions.condition_classification==1,:); % Get's the aluminum condition times
            sample_before_after = obj.Conditions.artifact_times/1000*250000; %Transfers the times that should be cut from the data to the audio sampling rate
                 % Find all times of aluminum event in data concatenate them
                % and create a vector for the samples to fit the sample in
                % the original file
            %data_aluminum = []; % an empty vector for data from the aluminum instances only
            samples_aluminum = []; % an empty vector for the sample times of the aluminum condition data
            for al = 1:length(alum_times) % iterates over the aluminum times to build two vectors: samples, data
                time_start = alum_times(al,1) + sample_before_after(2);
                time_end = alum_times(al,2) - sample_before_after(1);
                temp_sample = time_start:time_end;
                samples_aluminum = [samples_aluminum,temp_sample];  
            end
            %data_aluminum = data(samples_aluminum);
            %clear data
            %data_aluminum = wdenoise(double(data_aluminum),4);
            %data_aluminum = normalize(data_aluminum,'range',[-1,1]);
            
        end
        function [cams] = syncCamera(obj,x)
            %Sync the camera frame times to the open ephys strobe times
            %   Choose the camera file for face camera or for whisker camera and
            %   extract the frames, find the starting frame from both the camera csv
            %   file and the strobe data from the open ephys. 
            %   x - categorical. recieves 1 or 2, 1 for whisker camera, 2 for face cam
            % Choose if whisking camera or face camera and choose the appropriate file
            cd(obj.Pathways.main)
            switch x
                case 1
                    disp('Choose csv file for the whisking camera');
                    strobe_channel = obj.digital_channels.whisker_cam;
                case 2
                    disp('Choose csv file for the face camera');
                    strobe_channel = obj.digital_channels.face_cam;
            end
            [file,path] = uigetfile('*.csv');
            cam_fps = input(['Please input the fps for the current cam' '\n']);


            % read csv file to variable
            csv_data = csvread([path,'/',file]);
            % read openephys all timestamps and digital channels timestamps npy files
            full_timestamps = readNPY([obj.Pathways.exp,'/experiment1/recording1/continuous/',obj.Pathways.FPGA,'/timestamps.npy']);
            open_ephys_TTL = load_open_ephys_binary([obj.Pathways.exp,'/experiment1/recording1/structure.oebin'],'events',1);
            % find camera strobe on times from timestamps
            strobe_data = double(open_ephys_TTL.Timestamps(open_ephys_TTL.Data == strobe_channel));
            % find the frames from both openephys and csv after the forced stop
            strobe_data(2:end,2) =  (strobe_data(2:end,1) - strobe_data(1:end-1))/open_ephys_TTL.Header.sample_rate;
            csv_data(2:end,3) = (csv_data(2:end,2) - csv_data(1:end-1,2))/10000000;
            strobe_start = find(strobe_data(1:3*cam_fps,2) > (5*(1/cam_fps)));
            strobe_start = strobe_data(strobe_start(end),1);
            csv_start = find(csv_data(1:3*cam_fps,3) > (5*(1/cam_fps)));
            csv_start = csv_start(end);
            % align both time stamps
            csv_times = (csv_data(:,2) - csv_data(1,2))/10000000;
            csv_times = csv_times - csv_times(csv_start);
            csv_frames = floor(csv_times * open_ephys_TTL.Header.sample_rate);
            csv_frames = csv_frames + (strobe_start - double(full_timestamps(1)));
            % save the timestamps and frame times
            cams.strobe_start_frame = strobe_start - double(full_timestamps(1));
            cams.csv_start_frame = csv_start;
            cams.csv_aligned_frames = csv_frames;
            cams.fps = cam_fps;
        end
        
        function [evokedEvents] = getEvoked(obj)
           data =  load_open_ephys_binary([obj.Pathways.exp,'/experiment1/recording1/structure.oebin'],'continuous',1,'mmap');
           spk_data = data.Data.Data.mapped(32+obj.analog_channels.speaker,:);
           [b,a] = butter(4,[100]/30000,"high");
           spk_filt = filtfilt(b,a,double(spk_data(1:12000000)));
           clear data
           
           timesEv = detectEvents(spk_filt);
           %complexSound = input('Input 1 if the expereriment has a complex sound stimulation with frequancies and imitated whisking or 2 for simple sine');
           if length(timesEv) < 600
               maxFreq = 150;
           else
               maxFreq = 600;
           end
           try
               %timesEv(1) = [];
               freqVector = [1:15:maxFreq];
               freqLines = [0:14]';
               evokedEvents.eventsFreqKhz = [4,5.5,7.5,9.5,11,13,15,16.5,18.5,20,22,24,25.5,27.5,30]';
               freqMatrix = freqVector+freqLines;
               evokedEvents.whiskTimes = timesEv((maxFreq+1):end);
               evokedEvents.timePerFreq = timesEv(freqMatrix);
           catch
               'simple evoked stim'
           end
           evokedEvents.Alltimes = timesEv;
           evokedEvents.pVal = zeros(length(obj.Units.good.id)+length(obj.Units.mua.id),30);
           for fr=1:15
               evokedEvents.pVal(:,(fr*2)-1:(fr*2)) = raster4freq(obj,evokedEvents.timePerFreq(fr,:));
           end
           evokedEvents.pVal(:,(fr+1)*2-1:(fr+1)*2) = raster4freq(obj,evokedEvents.whiskTimes);
           for i =1:size(evokedEvents.pVal,1)
               intSum = sum(evokedEvents.pVal(i,2:2:32));
               evokedEvents.pVal(i,33) = (intSum > 0);
           end
           
           
           
           
            function timesEv = detectEvents(data)
                data = abs(data); 
                thr = 100;
                timesEv = find(data > thr);
%                 below_thr = find(data < thr)+1;
%                 timesEv = intersect(above_thr',below_thr');

                timeBefore = 5*30000;
                timesEv((timesEv - timeBefore < 1)) = [];
                timesEv([false,diff(timesEv) < 600]) = [];
%                 wind_mat = (ones(length(timesEv),timeBefore) .* [-timeBefore:1:-1]) + timesEv;
%                 wind_mat = round(wind_mat);
%                 wind_dat_onset = data(wind_mat);
%                 [rows,~] = find(wind_dat_onset > (thr)); rows = unique(rows);
%                 timesEv(rows) = [];

            end
            function sign4freq = raster4freq(obj,eventVec)
                for goodRast=1:length(obj.Units.good.id)
                    tempraster = createRaster(eventVec,obj.Units.good.times(goodRast,:),0.5);
                    midRast = length(tempraster)/2;
                    before = mean(tempraster(:,1:midRast),2);
                    after = mean(tempraster(:,midRast:midRast+3000),2);
                    [sign4freq(goodRast,1),sign4freq(goodRast,2)] = signrank(after,before,'tail','right');
                end
                if isempty(goodRast)
                    goodRast = 0;
                end
                for muaRast=1:length(obj.Units.mua.id)
                    try
                        midRast = length(tempraster)/2;
                        tempraster = createRaster(eventVec,obj.Units.mua.times(muaRast,:),0.5);
                    catch
                        tempraster = createRaster(eventVec,obj.Units.mua.times(muaRast,:),0.5);
                        midRast = length(tempraster)/2;
                    end
                    before = mean(tempraster(:,1:midRast),2);
                    after = mean(tempraster(:,midRast:midRast+3000),2);
                    [sign4freq(goodRast+muaRast,1),sign4freq(goodRast+muaRast,2)] = signrank(after,before,'tail','right');
                end                
            end
        end
        function conditionV = getConVector(obj)
%             artifactStart = obj.Conditions.artifact_times(2) * 30;
%             artifactEnd = obj.Conditions.artifact_times(1) * 30;
%             all_times = [obj.Conditions.all_changes + artifactStart,[obj.Conditions.all_changes(2:end);length(obj.sound_events.fullSound)]-artifactEnd];
            fullVector = mat2vec(obj.Conditions.all_changesClean);
            aluminumTimes = obj.Conditions.all_changesClean(obj.Conditions.classificationCleaned == 1,:);
            mutedTimes = obj.Conditions.all_changesClean(obj.Conditions.classificationCleaned == 2,:);
            nonTimes = obj.Conditions.all_changesClean(obj.Conditions.classificationCleaned == 3,:);
            % create idx vector for aluminum
            conditionV.aluminum = mat2vec(aluminumTimes);
            conditionV.muted = mat2vec(mutedTimes);
            conditionV.non = mat2vec(nonTimes);
            allWhisking = mat2vec(obj.Whisking.all);
            longWhisking = mat2vec(obj.Whisking.long);
            conditionV.noWhisking = fullVector(~ismember(fullVector,allWhisking));
            

            conditionV.aluminum = conditionV.aluminum(ismember(conditionV.aluminum,allWhisking));
            conditionV.non = conditionV.non(ismember(conditionV.non,allWhisking));
            conditionV.muted = conditionV.muted(ismember(conditionV.muted,allWhisking));
            conditionV.fullVector = fullVector;
            
            
            function collapsedVector = mat2vec(startEndMat)
                collapsedVector = [];
                for epoch=1:length(startEndMat)
                    tempVec = [startEndMat(epoch,1):startEndMat(epoch,2)];
                    collapsedVector = [collapsedVector,tempVec];
                end
            end

              
            
        end

        function [unitTimesCon] = unitTimesCondition(obj)
            
            alum_temp_good = getPerCon(obj.Units.good.times,obj.conditionVector.aluminum);
            alum_temp_mua = getPerCon(obj.Units.mua.times,obj.conditionVector.aluminum);
            
            unitTimesCon.aluminum.times = [alum_temp_good.times;alum_temp_mua.times];
            unitTimesCon.aluminum.FR = [alum_temp_good.FR;alum_temp_mua.FR];
            
            
            mut_temp_good = getPerCon(obj.Units.good.times,obj.conditionVector.muted);
            mut_temp_mua = getPerCon(obj.Units.mua.times,obj.conditionVector.muted);
            
            unitTimesCon.muted.times = [mut_temp_good.times;mut_temp_mua.times];
            unitTimesCon.muted.FR = [mut_temp_good.FR;mut_temp_mua.FR];
            
            non_temp_good = getPerCon(obj.Units.good.times,obj.conditionVector.non);
            non_temp_mua = getPerCon(obj.Units.mua.times,obj.conditionVector.non);
            
            unitTimesCon.non.times = [non_temp_good.times;non_temp_mua.times];
            unitTimesCon.non.FR = [non_temp_good.FR;non_temp_mua.FR];
            
            noWhisking_temp_good = getPerCon(obj.Units.good.times,obj.conditionVector.noWhisking);
            noWhisking_temp_mua = getPerCon(obj.Units.mua.times,obj.conditionVector.noWhisking);
            
            unitTimesCon.noWhisking.times = [noWhisking_temp_good.times;noWhisking_temp_mua.times];
            unitTimesCon.noWhisking.FR = [noWhisking_temp_good.FR;noWhisking_temp_mua.FR];
            
            
            
            
            
            function condMat = getPerCon(units,conVector)
                condMat.times = cell(size(units,1),1);
                isCon = ismember(units,conVector);
                condMat.FR = sum(isCon,2)/(length(conVector)/30000);
                for unit = 1:size(units)
                    condMat.times{unit} = units(unit,isCon(unit,:));
                end
            end
        end
        
        
        function conditionStats = getBinnedStats(obj,binSize)
            unitsCon = obj.unitTimesCondition;
            binSizeSF = binSize * obj.SF;
            alumMatrix = allBinMat(obj.conditionVector.aluminum',binSizeSF);
            alumBinFR = getBinnedFR(unitsCon.aluminum.times,alumMatrix)/binSize;
            mutedMatrix = allBinMat(obj.conditionVector.muted',binSizeSF);
            mutedBinFR = getBinnedFR(unitsCon.muted.times,mutedMatrix)/binSize;
            nonMatrix = allBinMat(obj.conditionVector.non',binSizeSF);
            nonBinFR = getBinnedFR(unitsCon.non.times,nonMatrix)/binSize;
            noWhiskMatrix = allBinMat(obj.conditionVector.noWhisking',binSizeSF);
            noWhiskBinFR = getBinnedFR(unitsCon.noWhisking.times,noWhiskMatrix);
            conditionStats.unitsPerCon = unitsCon;
            conditionStats.binned.aluminum = alumBinFR;
            conditionStats.binned.muted = mutedBinFR;
            conditionStats.binned.non = nonBinFR;
            conditionStats.binned.noWhisk = noWhiskBinFR;
            conditionStats.STD(:,1) = std(alumBinFR,0,2);
            conditionStats.STD(:,2) = std(mutedBinFR,0,2);
            conditionStats.STD(:,3) = std(nonBinFR,0,2);
            conditionStats.STD(:,4) = std(noWhiskBinFR,0,2);
            

            function binMatrixFull = allBinMat(condVec,binSize)
                condParsed = find(diff(condVec) > 1);
                condParsedAll = [1;condParsed+1];
                condParsedAll(:,2) = [condParsed;length(condVec)];
                fittedLength = zeros(size(condParsedAll,1),1);
                for epoch = 1:size(condParsedAll,1)
                    epochLength = (condParsedAll(epoch,2) - condParsedAll(epoch,1) + 1);
                    fittedLength(epoch) = floor(epochLength/binSize)*binSize;
                    condParsedAll(epoch,2) = fittedLength(epoch) + condParsedAll(epoch,1) - 1;
                end
                allBins = sum(fittedLength/binSize);
                binMatrixFull = zeros(allBins,binSize);
                k = 1;
                for epochBins = 1:size(condParsedAll,1)
                    binnedVec = [condParsedAll(epochBins,1):condParsedAll(epochBins,2)];
                    binMatrix = reshape(binnedVec,binSize,[])';
                    raster4freq(binMatrixFull(k:k+size(binMatrix,1)-1,:)) = binMatrix;
                    k = size(binMatrix,1)+k;
                    
                end 
                binMatrixFull = condVec(binMatrixFull);
            end
            
            function binnedFR = getBinnedFR(conUnits,binnedMat)
               binnedFR = zeros(size(binnedMat,1),length(conUnits));
               for unit = 1:length(conUnits) 
                   tempFR = ismember(binnedMat,conUnits{unit});
                   binnedFR(:,unit) = sum(tempFR,2);
                   
               end
               binnedFR = binnedFR';
            end
            
        end

        function [SDFstruct] = getSDF(obj)
            expLength = length(obj.sound.full);
            expLengthDS = expLength/30000;
            SDFstruct.timeDS = [1:30:expLength]';
            conditions = fieldnames(obj.conditionVector);
            for cn = 1:length(conditions)
                SDFstruct.vecDS.(conditions{cn}) = intersect(SDFstruct.timeDS,obj.conditionVector.(conditions{cn}));
            end
            SDFstruct.expLength = expLength;
            SDFstruct.expLengthDS = expLengthDS;
            goodUnits = length(obj.Units.good.id);
            muaUnits = length(obj.Units.mua.id);
            unitsNum = goodUnits + muaUnits;
            SDFstruct.units = ones(unitsNum,length(SDFstruct.timeDS));
            SDFstruct.raster = ones(unitsNum,length(SDFstruct.timeDS));
            for goodU = 1:goodUnits
                tempUnit = nonzeros(obj.Units.good.times(goodU,:))/30000;
                [SDFstruct.units(goodU,:),SDFstruct.raster(goodU,:)] = getSDFBinned(tempUnit,expLengthDS);
            end
            for muaU = 1:muaUnits
                tempUnit = nonzeros(obj.Units.mua.times(muaU,:))/30000;
                [SDFstruct.units(goodU+muaU,:),SDFstruct.raster(goodU+muaU,:)] = getSDFBinned(tempUnit,expLengthDS);
            end
            SDFstruct.soundSmoothed = resample(double(obj.sound.smoothed),1000,30000);
        end

        
        
        
  
        

        
        
    end
end

