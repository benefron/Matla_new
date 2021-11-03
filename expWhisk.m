classdef expWhisk
    %A class That contain the information and the methods for analyzing the
    %auditory whisking experiments
    %   Detailed explanation goes here
    
    properties
        expName 
        SF = 30000; % defult SF can be changed if needed
        whisking_times
        touch_times
        sound_times
        expType % 1 - for awake 2- for anestheized galvo
        clusterSpikes
        meanSound
        meanWhisk
        meanObjectTouch
        windowSize = 2; % defult window size for all analysis
        expLength 
        statistics
        statWind = 0.5; % default window for calculating statistics
        filterCoF
        conditionsTimes
        
    end

    properties(Hidden = 1)
        openEphysDataType
        soundChannel
        whiskingChannel
        objectTouchChannel
        spkChannel
        paramsFile
        exCD
        color
        
    end
        
    
    methods
        function obj = expWhisk(params,exName,varargin)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.expName = exName;
            if nargin > 2
                obj.windowSize = varargin{1};
            end
            funCase = 4;
            rerun = 1;
            run(params);
            obj.SF = SF;
            if openEphDataType == 3
                obj.conditionsTimes.condition_1 = [expTimes(1,1),expTimes(1,2)];
                obj.conditionsTimes.condition_2 = [expTimes(2,1),expTimes(2,2)];
                obj.conditionsTimes.condition_3 = [expTimes(3,1),expTimes(3,2)];
            elseif openEphDataType == 2
                obj.conditionsTimes = extRandConditions(params,binData);
            else
                obj.conditionsTimes = extRandConditions(params);
            end
            [whisking_times] = Whisking_times_per_condition(params,1);
            %[sound_times] = Whisking_times_per_condition(params,6);
            if openEphDataType == 2
                touch_time = touch_times(params,binData);
            else
                touch_time = touch_times(params);
            end
            touch_time.noObject = whisking_times.noObject;
            if evokedRes == 1
                touch_time.Evoked = whisking_times.Evoked;
            end
            obj.expType = evokedRes;
            obj.clusterSpikes = sorted_from_phy(params);
            obj.expLength = max([max(max(obj.clusterSpikes.good.times)),max(max(obj.clusterSpikes.mua.times))]);
            obj.whisking_times = obj.windowed(whisking_times);
            obj.touch_times = obj.windowed(touch_time);
            %obj.sound_times = obj.windowed(sound_times);
            obj.openEphysDataType = openEphDataType;
            obj.soundChannel = soundChannel;
            obj.whiskingChannel = whiskingChannel;
            obj.objectTouchChannel = objectTouchChannel;
            obj.spkChannel = speakerChannel;
            obj.paramsFile = params;
            obj.exCD = exCD;
            obj.filterCoF.b = b;
            obj.filterCoF.a = a;
            if openEphDataType == 2
                obj.meanWhisk = obj.getMeans(1,binData);
                obj.meanObjectTouch = obj.getMeans(2,binData);
                obj.meanSound = obj.getMeans(3,binData);
            else
                obj.meanWhisk = obj.getMeans(1);
                obj.meanObjectTouch = obj.getMeans(2);
                obj.meanSound = obj.getMeans(3);
            end
            obj.statistics = obj.calc_statistics();
            obj.color.Aluminum = [0.0039 * 210 , 0.0039 * 30 , 0.0039 * 75]; % D21E4B adobe color wheel
            obj.color.noObject = [0.0039 * 232 , 0.0039 * 134 , 0.0039 * 65]; % 284B59 adobe color wheel
            obj.color.Muted = [0.0039 * 40, 0.0039 * 75, 0.0039 * 89]; % E88641 adobe color wheel
            
        end
        
        function meanArt = getMeans(obj,meanADC,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            cd(obj.exCD)
            switch meanADC
                case 1 % whisk
                    switch obj.openEphysDataType
                        case 1 % openephys format
                        dataRaw = load_open_ephys_data_faster(obj.whiskingChannel)';
                        case 2 % .dat bin file format
                        dataRaw = varargin{1}.mapped(obj.whiskingChannel,1:end);
                        case 3
                        dataRaw = obj.whiskingChannel';
                    end
                    % Smoothing the Data for the detection
                    [b,a] = butter(2,[20,200]/(obj.SF/2),'bandpass');
                    data = double(dataRaw);
                    data = filtfilt(b,a,data);
                    data = normalize(data);
                    %sigmaTimeSec = 0.01;
                    %sigmaSamples = round(sigmaTimeSec*obj.SF);
                    %window = fspecial('gaussian',[1, sigmaSamples*6],sigmaSamples);
                    %data = conv(data,window,'same');
                     
                case 2 %objectTouch
                    switch obj.openEphysDataType
                        case 1 % openephys format
                        dataRaw = load_open_ephys_data_faster(obj.objectTouchChannel)';
                        case 2 % .dat bin file format
                        dataRaw = varargin{1}.mapped(obj.objectTouchChannel,1:end);
                        case 3
                        dataRaw = obj.objectTouchChannel';
                    end
                    
                    dataAll = double(dataRaw);
                    dataAll = filtfilt(obj.filterCoF.b,obj.filterCoF.a,dataAll);
                    dataAlum = [];
                    dataMut = [];
                    for vec = 1:size(obj.conditionsTimes.condition_1,1)
                        tempCond = [obj.conditionsTimes.condition_1(vec,1):obj.conditionsTimes.condition_1(vec,2)];
                        dataAlum = [dataAlum,tempCond];
                    end
                    
                    for vecM = 1:size(obj.conditionsTimes.condition_2,1)
                        tempCond = [obj.conditionsTimes.condition_2(vecM,1):obj.conditionsTimes.condition_2(vecM,2)];
                        dataMut = [dataMut,tempCond];
                    end
                    
                    dataAlum(2,:) = dataAll(dataAlum(1,:));
                    dataMut(2,:) = dataAll(dataMut(1,:));
                    dataAlum(2,:) = normalize(dataAlum(2,:));
                    dataMut(2,:) = normalize(dataMut(2,:));
                    data = dataAll;
                    data(dataAlum(1,:)) = dataAlum(2,:);
                    data(dataMut(1,:)) = dataMut(2,:);
                    
                    
                     
                case 3
                    data = audioread(obj.soundChannel);
                    if obj.expType == 1
                        switch obj.openEphysDataType
                            case 1 % openephys format
                                dataSPK = load_open_ephys_data_faster(obj.spkChannel)';
                            case 2 % .dat bin file format
                                openEphDataType = 2; funCase = 99; rerun = 1;
                                run(obj.paramsFile)
                                dataSPK = binData.mapped(obj.spkChannel,1:end);
                            case 3
                                dataSPK = obj.spkChannel';
                        end
                    end
                    data = abs(data);
                    
                    
            end
            fn = fieldnames(obj.whisking_times);
            for con = 1:length(fn)
                meanArt.(fn{con}).whisk = mean(data(obj.whisking_times.(fn{con})(:,2:end-2)));
                if con < 3
                    meanArt.(fn{con}).touch = mean(data(obj.touch_times.(fn{con})(:,2:end-2)));
                end
                %meanArt.(fn{con}).sound = mean(data(obj.sound_times.(fn{con})));
            end
            if meanADC == 3
                if length(fn) > 3
                    con = 4;
                    meanArt.(fn{con}) = mean(dataSPK(obj.whisking_times.(fn{con})(:,5:end-5)));
                end
            end
            
        end
        
        function windowedWhisk = windowed(obj,whisking_times)
            % A function to create the windoed times around the whisking
            % times
            windowVector = (-obj.windowSize/2*obj.SF+1):1:(obj.windowSize/2*obj.SF);
            fn = fieldnames(whisking_times);
            for con = 1:length(fn)
                whisking_times.(fn{con})((whisking_times.(fn{con}) - obj.windowSize * obj.SF <= 0 ) | (whisking_times.(fn{con}) + obj.windowSize * obj.SF > obj.expLength )) = [];
                windowedWhisk.(fn{con}) = (ones(length(whisking_times.(fn{con})),length(windowVector))).* windowVector + whisking_times.(fn{con});
            end
            
            
        end
        
        function stats = calc_statistics(obj)
            clusterID = [obj.clusterSpikes.good.id;obj.clusterSpikes.mua.id];
            classification = [char(ones(length(obj.clusterSpikes.good.id),1) * 'good');char(ones(length(obj.clusterSpikes.mua.id),1) * 'mua ')];
            if obj.expType == 1
                evoked = signCalc(obj,'Evoked');
            else
                evoked = ones(length(clusterID),2);
            end
            A_A = signCalc(obj,'Aluminum');
            M_M = signCalc(obj,'Muted');
            N_N = signCalc(obj,'noObject');
            A_M = rankCalc(obj,{'Aluminum';'Muted'});
            A_N = rankCalc(obj,{'Aluminum';'noObject'});
            M_N = rankCalc(obj,{'Muted';'noObject'});
            Experimant = char(ones(length(clusterID),1) * obj.expName);
            stats = table(Experimant,clusterID,classification,evoked,A_A,M_M,N_N,A_M,A_N,M_N);
            
            function colName = signCalc(obj,condition)
                zeroTime = 0.5 * obj.windowSize * obj.SF;
                statWindTemp = obj.statWind * obj.SF;
                startCalc = zeroTime - 1 - statWindTemp;
                endCalc = zeroTime + 0.2 * obj.SF;
                colName = zeros(length(clusterID),2);
                timeAfter = endCalc - zeroTime;
                ind = 1;
                for good =1:length(obj.clusterSpikes.good.id)
                    cluster = ismember(obj.touch_times.(condition),obj.clusterSpikes.good.times(good,:));
                    [colName(ind,1),colName(ind,2)] = signrank(sum(cluster(:,startCalc:zeroTime - 1),2)./statWindTemp,sum(cluster(:,zeroTime:endCalc),2)./timeAfter);
                    ind = ind + 1;
                end
                for mua = 1:length(obj.clusterSpikes.mua.id)
                    cluster = ismember(obj.touch_times.(condition),obj.clusterSpikes.mua.times(mua,:));
                    [colName(ind,1),colName(ind,2)] = signrank(sum(cluster(:,startCalc:zeroTime - 1),2)./statWindTemp,sum(cluster(:,zeroTime:endCalc),2)./timeAfter);
                   ind = ind + 1;
                end
            end
            function colName = rankCalc(obj,conditions)
                zeroTime = 0.5 * obj.windowSize * obj.SF;
                statWindTemp = obj.statWind * obj.SF;
                startCalc = zeroTime - 1 - statWindTemp;
                endCalc = zeroTime + 0.2 * obj.SF;
                colName = zeros(length(clusterID),2);
                timeAfter = endCalc - zeroTime;
                ind = 1;
                for good =1:length(obj.clusterSpikes.good.id)
                    cond1 = ismember(obj.touch_times.(conditions{1}),obj.clusterSpikes.good.times(good,:));
                    cond2 = ismember(obj.touch_times.(conditions{2}),obj.clusterSpikes.good.times(good,:));
                    cond1 = sum(cond1(:,zeroTime:endCalc),2)./timeAfter - sum(cond1(:,startCalc:zeroTime - 1),2)./statWindTemp;
                    cond2 = sum(cond2(:,zeroTime:endCalc),2)./timeAfter - sum(cond2(:,startCalc:zeroTime - 1),2)./statWindTemp;
                    [colName(ind,1),colName(ind,2)] = ranksum(cond1,cond2);
                    ind = ind + 1;
                end
                for mua = 1:length(obj.clusterSpikes.mua.id)
                    cond1 = ismember(obj.touch_times.(conditions{1}),obj.clusterSpikes.mua.times(mua,:));
                    cond2 = ismember(obj.touch_times.(conditions{2}),obj.clusterSpikes.mua.times(mua,:));
                    cond1 = sum(cond1(:,zeroTime:endCalc),2)./timeAfter - sum(cond1(:,startCalc:zeroTime - 1),2)./statWindTemp;
                    cond2 = sum(cond2(:,zeroTime:endCalc),2)./timeAfter - sum(cond2(:,startCalc:zeroTime - 1),2)./statWindTemp;
                    [colName(ind,1),colName(ind,2)] = ranksum(cond1,cond2);
                    ind = ind + 1;
                end
            end
        end
        
        function clustFig = createFigure(obj,conditions,cluster,psBin,varargin)
            stimType = 'touch_times';
            stimNum = 1;
            
            if ~isempty(varargin)
                stimType = 'whisking_times';
            end
            if conditions{1}(1) == 'A'
                if conditions{2}(1) == 'M'
                    statCol = 'A_M';
                    figureCond = 1;
                else
                    statCol = 'A_N';
                    figureCond = 2;
                end
            else
                statCol = 'M_N';
            end
            class = obj.statistics.classification(find(obj.statistics.clusterID == cluster),:);
            
            if size(class,1) > 1
                inClass = input('input the cluster classification: 1 for "good", 2 for "mua" ')
                class = class(inClass,:);
            end
            if class == 'mua '
                class = class(1:3);
            end
            clusterInd = find(obj.clusterSpikes.(class).id == cluster);
            tableInd = find(obj.statistics.clusterID == cluster);
            if length(tableInd) > 1
                tableInd = tableInd(inClass);
            end
            switch figureCond
                case 1
                    cond1_cluster =  ismember(obj.(stimType).(conditions{1}),obj.clusterSpikes.(class).times(clusterInd,:));
                    cond2_cluster =  ismember(obj.(stimType).(conditions{2}),obj.clusterSpikes.(class).times(clusterInd,:));
                case 2
                    cond1_cluster =  ismember(obj.(stimType).(conditions{1}),obj.clusterSpikes.(class).times(clusterInd,:));
                    cond2_cluster =  ismember(obj.whisking_times.(conditions{2}),obj.clusterSpikes.(class).times(clusterInd,:));            
            end
            clustFig = figure('units','normalized');
            ax1 = subplot(3,1,1);
            getRaster(obj,cond1_cluster,obj.color.(conditions{1}));
            set(gca,'xtick',[])
            title(conditions{1})
            ylim([1 200]);
            ax2 = subplot(3,1,2);
            getRaster(obj,cond2_cluster,obj.color.(conditions{2}));
            set(gca,'xtick',[])
            ylim([1 200]);
            title(conditions{2})
            ax3 = subplot(3,1,3);
            [time1,vec1] = Psth_vec(sum(cond1_cluster), psBin , obj.SF );
            normFactor1 = size(cond1_cluster,1);
            [time2,vec2] = Psth_vec(sum(cond2_cluster), psBin , obj.SF );
            normFactor2 = size(cond2_cluster,1);
            plot(time1/1000-0.5 * obj.windowSize,smooth(vec1/normFactor1,'sgolay',3),'Color',obj.color.(conditions{1}));
            hold on
            b2 = plot(time2/1000-0.5 * obj.windowSize,smooth(vec2/normFactor2,'sgolay',3),'Color',obj.color.(conditions{2}));
            %b2.FaceAlpha = 0.85;
            %set(gca,'xtick',[])
            t = title(['p = ' num2str(obj.statistics.(statCol)(tableInd,1))]);
%             ax4 = subplot(4,1,4);
%             if ~isempty(varargin)
%                 p1 = plot((1:length(obj.meanWhisk.(conditions{1}).whisk))/obj.SF - obj.windowSize/2,obj.meanWhisk.(conditions{1}).touch-mean(obj.meanWhisk.(conditions{1}).whisk),'color',obj.color.(conditions{1}));
%                 hold on
%                 p2 = plot((1:length(obj.meanWhisk.(conditions{2}).whisk))/obj.SF - obj.windowSize/2,obj.meanWhisk.(conditions{2}).touch-mean(obj.meanWhisk.(conditions{2}).whisk),'color',obj.color.(conditions{2}));
%                 p2.Color(4) = 0.9;
%             else
%                 p1 = plot((1:length(obj.meanObjectTouch.(conditions{1}).touch))/obj.SF - obj.windowSize/2,obj.meanObjectTouch.(conditions{1}).touch-mean(obj.meanObjectTouch.(conditions{1}).touch),'color',obj.color.(conditions{1}));
%                 hold on
%                 p2 = plot((1:length(obj.meanObjectTouch.(conditions{2}).touch))/obj.SF - obj.windowSize/2,obj.meanObjectTouch.(conditions{2}).touch-mean(obj.meanObjectTouch.(conditions{2}).touch),'color',obj.color.(conditions{2}));
%                 p2.Color(4) = 0.9;
%             end
            %set(gca,'xtick',[])
            title('ObjectTouch')
%             ax5 = subplot(6,1,5);
%             if obj.expType == 0
%                  t1 = plot((1:length(obj.meanWhisk.(conditions{1})))/obj.SF - obj.windowSize/2,obj.meanWhisk.(conditions{1}),'color',obj.color.(conditions{1}));
%                  hold on
%                  t2 = plot((1:length(obj.meanWhisk.(conditions{2})))/obj.SF - obj.windowSize/2,obj.meanWhisk.(conditions{2}),'color',obj.color.(conditions{2}));
%                  title('Galvo signal')
%                  t2.Color(4) = 0.8;
%             else
%                  t1 = plot((1:length(obj.meanObjectTouch.(conditions{1}).touch))/obj.SF - obj.windowSize/2,obj.meanWhisk.(conditions{1}).touch,'color',obj.color.(conditions{1}));
%                  hold on
%                  if conditions{2}(1) == 'n'
%                      t2 = plot((1:length(obj.meanObjectTouch.(conditions{2}).whisk))/obj.SF - obj.windowSize/2,obj.meanWhisk.(conditions{2}).whisk,'color',obj.color.(conditions{2}));
%                  else
%                      t2 = plot((1:length(obj.meanObjectTouch.(conditions{2}).touch))/obj.SF - obj.windowSize/2,obj.meanWhisk.(conditions{2}).touch,'color',obj.color.(conditions{2}));
%                  end
%                  title('Whisking')
%                  t2.Color(4) = 0.8;
%             end
%             ax6 = subplot(6,1,6)
%             so1 = plot((1:length(obj.meanObjectTouch.(conditions{1}).touch))/obj.SF - obj.windowSize/2,obj.meanSound.(conditions{1}).touch-mean(obj.meanSound.(conditions{1}).touch),'color',obj.color.(conditions{1}));
%             hold on
%             if figureCond == 1
%                 so2 = plot((1:length(obj.meanObjectTouch.(conditions{2}).touch))/obj.SF - obj.windowSize/2,obj.meanSound.(conditions{2}).touch-mean(obj.meanSound.(conditions{2}).touch),'color',obj.color.(conditions{2}));
%             else
%                 so2 = plot((1:length(obj.meanSound.(conditions{2}).whisk))/obj.SF - obj.windowSize/2,obj.meanSound.(conditions{2}).whisk-mean(obj.meanSound.(conditions{2}).whisk),'color',obj.color.(conditions{2}));
% 
%             end
%             so2.Color(4) = 0.3;
%             title('Sound')
%             
            sgt = sgtitle({
                [obj.expName]
                ['Cluster #' num2str(cluster)]
                [conditions{1} ' vs ' conditions{2}]
                });
            sgt.FontSize = 15;
            linkaxes([ax1,ax2,ax3],'x')
            %xlim([-0.3 0.8])
       
            
            
            
            
            
            function getRaster(obj,condition,rastColor)
                hold on
                for ind = 1:size(condition,1)
                    y = ind*ones(1,length(find(condition(ind,:) == 1)));
                    scatter((find(condition(ind,:) == 1)/obj.SF)- 0.5 * obj.windowSize,y,[],rastColor,'.')
                end
            end
            
        end
        
        function evoFig = evokedFigure(obj,cluster)
            class = obj.statistics.classification(find(obj.statistics.clusterID == cluster),:);
            if size(class,1) > 1
                inClass = input('input the cluster classification: 1 for "good", 2 for "mua" ')
                class = class(inClass,:);
            end
            if class == 'mua '
                class = class(1:3)
            end
            clusterInd = find(obj.clusterSpikes.(class).id == cluster);
            cond_cluster =  ismember(obj.whisking_times.Evoked,obj.clusterSpikes.(class).times(clusterInd,:));
            evoFig = figure; 
            ax1 = subplot(2,1,1);
            getRaster(obj,cond_cluster,'k');
            ax2 = subplot(2,1,2);
            plot((((1:length(obj.meanSound.Evoked))/obj.SF)-0.5 * obj.windowSize),obj.meanSound.Evoked,'k')
            
            function getRaster(obj,condition,rastColor)
                hold on
                for ind = 1:size(condition,1)
                    y = ind*ones(1,length(find(condition(ind,:) == 1)));
                    scatter((find(condition(ind,:) == 1)/obj.SF)- 0.5 * obj.windowSize,y,[],rastColor,'.')
                end
            end
            linkaxes([ax1,ax2],'x');
        end
        
        function STA = createSTA(obj,conditions,cluster,varargin)
            funCase = 99;
            rerun = 1;
            openEphDataType = obj.openEphysDataType;
            run(obj.paramsFile)
            class = obj.statistics.classification(find(obj.statistics.clusterID == cluster),:);
            if size(class,1) > 1
                inClass = input('input the cluster classification: 1 for "good", 2 for "mua" ')
                class = class(inClass,:);
            end
            if class == 'mua '
                class = class(1:3);
            end
            clusterInd = find(obj.clusterSpikes.(class).id == cluster);
            whiskStart = obj.windowSize * obj.SF/2 - 0.1 * obj.SF;
            whiskEnd = whiskStart + obj.windowSize * obj.SF/4;
            spikeTimes_1 = ismember(obj.clusterSpikes.(class).times(clusterInd,:),obj.whisking_times.(conditions{1})(:,whiskStart:whiskEnd));
            tempSP = obj.clusterSpikes.(class).times(clusterInd,:);
            spikeTimes_2 = ismember(obj.clusterSpikes.(class).times(clusterInd,:),obj.whisking_times.(conditions{2})(:,whiskStart:whiskEnd));
            spikeTimes_1 = tempSP(spikeTimes_1);
            spikeTimes_2 = tempSP(spikeTimes_2);
            cd(obj.exCD)
            data = audioread(obj.soundChannel);
            
            switch obj.openEphysDataType
                        case 1 % openephys format
                        dataRaw = load_open_ephys_data_faster(obj.whiskingChannel)';
                        case 2 % .dat bin file format
                        dataRaw = binData.mapped(obj.whiskingChannel,1:end);
                        case 3
                        dataRaw = obj.whiskingChannel';
            end
                    % Smoothing the Data for the detection
                    [b,a] = butter(4,2/(obj.SF/2),'high');
                    dataWhisk = double(dataRaw);
                    dataWhisk = filtfilt(b,a,dataWhisk);
                    sigmaTimeSec = 0.01;
                    sigmaSamples = round(sigmaTimeSec*obj.SF);
                    window = fspecial('gaussian',[1, sigmaSamples*6],sigmaSamples);
                    dataWhisk = conv(dataWhisk,window,'same');
            
            
            windMat_1 = zeros(length(spikeTimes_1),obj.windowSize * obj.SF) + (1:obj.windowSize * obj.SF) - obj.windowSize * obj.SF/2 + spikeTimes_1';
            windMat_2 = zeros(length(spikeTimes_2),obj.windowSize * obj.SF) + (1:obj.windowSize * obj.SF) - obj.windowSize * obj.SF/2 + spikeTimes_2';
            mat1 = repmat((1:obj.windowSize * obj.SF)/obj.SF - 1,size(windMat_1,1),1);
            mat2 = repmat((1:obj.windowSize * obj.SF)/obj.SF - 1,size(windMat_2,1),1);
            w1 = sort(mat1(ismember(windMat_1,tempSP)));
            w2 = sort(mat2(ismember(windMat_2,tempSP)));
            n1 = size(mat1,1);
            n2 = size(mat2,1);
            
            
            soundVec = envelope(data-mean(data));
            figure; 
            p1 = subplot(3,1,1);pp1 = plot((1:obj.windowSize * obj.SF)/obj.SF - 1,mean(soundVec(windMat_1))-mean(mean(soundVec(windMat_1))));
            title('Sound');
            yl = ylim;
            hold on
            pp1.Color = obj.color.(conditions{1});
            pp2 = plot((1:obj.windowSize * obj.SF)/obj.SF - 1,mean(soundVec(windMat_2))-mean(mean(soundVec(windMat_2))));
            ylim(yl);
            pp2.Color = obj.color.(conditions{2});
            pp2.Color(4) = 0.1;
            p3 = subplot(3,1,3);
            pd1 = fitdist(w1,'kernel');
            pd2 = fitdist(w2,'kernel');
            xfit = (0 - obj.windowSize/2):0.005:(0 + obj.windowSize/2);
            histogram(w1,500,'FaceColor',obj.color.(conditions{1}),'EdgeColor','none','FaceAlpha',1,'Normalization','probability')           
            hold on
            histogram(w2,500,'FaceColor',obj.color.(conditions{2}),'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability')
            title('Spikes')
            
            fit1 = plot(xfit,pdf(pd1,xfit));%/max(pdf(pd1,xfit)));
            fit2 = plot(xfit,pdf(pd2,xfit));%/max(pdf(pd2,xfit)));
            fit1.Color = obj.color.(conditions{1});
            fit2.Color = obj.color.(conditions{2});
            legend([conditions{1} '- ' num2str(n1)],[conditions{2} '- ' num2str(n2)],[conditions{1} ' KDE'],[conditions{1} 'KDE'])
            p2 = subplot(3,1,2);
            ppW1 = plot((1:obj.windowSize * obj.SF)/obj.SF - 1,mean(dataWhisk(windMat_1))-mean(mean(dataWhisk(windMat_1))));
            title('Whisking');
            hold on
            ppW1.Color = obj.color.(conditions{1});
            ppW2 = plot((1:obj.windowSize * obj.SF)/obj.SF - 1,mean(dataWhisk(windMat_2))-mean(mean(dataWhisk(windMat_2)))); 
            ppW2.Color = obj.color.(conditions{2});

            STA.(conditions{1}) = windMat_1;
            STA.(conditions{2}) = windMat_2;
            
            
            sgt = sgtitle({
                ['STA of sound']
                [obj.expName]
                ['Cluster #' num2str(cluster)]
                [conditions{1} ' vs ' conditions{2}]
                });
            sgt.FontSize = 11;
            linkaxes([p1 p2 p3],'x')
            xlim([-1 1])
            

            

            
        end
        
        
        function [psth_mat,order_vec] = populationFiring(obj,condition,psBin,stimType)
            class = {'good','mua'};
            for classInd = 1:2
                %if obj.expType == 1
                    clusIT = obj.statistics.clusterID((obj.statistics.A_M(:,2) == 1 ) & obj.statistics.classification == class{classInd}(:,1));
                    clusI = find(ismember(obj.clusterSpikes.(class{classInd}).id,clusIT));
                %else
                    %clusI = obj.clusterSpikes.(class{classInd}).id;
                %end
                for cluster = 1:length(clusI)
                    clustTemp = ismember(obj.(stimType).(condition),obj.clusterSpikes.(class{classInd}).times(clusI(cluster),:));
                    [timeP,vec] = Psth_vec(sum(clustTemp),psBin,obj.SF);
                    normFactor = size(clustTemp,1);
                    psth_mat.(class{classInd})(cluster,:) = vec/normFactor;
                    %orderVec(cluster,1) = (sum(sum(clustTemp(:,obj.windowSize*obj.SF/2:floor(obj.windowSize*obj.SF/2+0.05*obj.SF))))-sum(sum(clustTemp(:,1:0.05*obj.SF)))/normFactor);
                end
                %[~,orderInd] = sort(orderVec,'descend');
                %psth_mat.(class{classInd}) = psth_mat.(class{classInd})(orderInd,:);
                %figure;
                %clims = [-1,1];
                %imagesc(timeP/1000-0.5*obj.windowSize,[1:size(psth_mat.(class{classInd}),1)],psth_mat.(class{classInd}),clims);
                %title([obj.expName, '- ',condition, ' ' , class{classInd}])
                %colorbar
                if length(clusI) < 1
                    psth_mat.(class{classInd}) = [];
                    order_vec.(class{classInd}) = [];
                else
                    ind = length(psth_mat.(class{classInd}));
                    order_vec.(class{classInd}) = mean(psth_mat.(class{classInd})(:,ind/2:ind/2+5),2)-mean(psth_mat.(class{classInd})(1:2),2);
                    order_vec.(class{classInd})(:,2) = clusI;
                end
                if exist('timeP')
                    if ~isempty(timeP)
                        psth_mat.time = timeP;
                    end
                end
            end

            
            
        end
        
        function [a,b] = populationFig(obj,psBin,stimType)
            [psth_alm,orderVec_alm] = obj.populationFiring('Aluminum',psBin,stimType);
            [psth_mut,orderVec_mut] = obj.populationFiring('Muted',psBin,stimType);
            class = {'good','mua'};
            for classInd = 1:2
                [~,orderInd] = sort(orderVec_alm.(class{classInd}),'descend');%+orderVec_mut.(class{classInd}))/2,'descend');
                psth_alm.(class{classInd}) = psth_alm.(class{classInd})(orderInd,:);
                psth_mut.(class{classInd}) = psth_mut.(class{classInd})(orderInd,:);
                figure
                
                a = subplot(1,2,1);
                imagesc(psth_alm.time/1000-0.5*obj.windowSize,[1:size(psth_alm.(class{classInd}),1)],psth_alm.(class{classInd}));
                title('Aluminum')
                
                b = subplot(1,2,2);
                imagesc(psth_mut.time/1000-0.5*obj.windowSize,[1:size(psth_mut.(class{classInd}),1)],psth_mut.(class{classInd}));
                title('Muted')
                colorbar
                sgtitle([obj.expName, '- ', class{classInd}])

                
            end
            
        end
        
            
        
       

    end
end

