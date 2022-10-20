function [all_psth] = mean_psth(event_times,unit_struct,time_window,binSize)

good_unit_times = unit_struct.good.times;
mua_unit_times = unit_struct.mua.times;

all_psth.good = zeros(size(good_unit_times,1),time_window*1000/binSize);
all_psth.mua = zeros(size(mua_unit_times,1),time_window*1000/binSize);
bEvent = 0;
aEvent = 0.075;
baseLine = bEvent + aEvent;
binNum4Stat = baseLine*time_window*30000;
all_psth.stats = [];
eventPoint = time_window/2*30000;

beforeEvent = eventPoint - bEvent*30000;
afterEvent = eventPoint + aEvent*30000 - 1;


for good = 1:size(good_unit_times,1)
   temp_raster = createRaster(event_times,good_unit_times(good,:),time_window);
   [all_psth.stats(good,1),all_psth.stats(good,2)] = signrank(sum(temp_raster(:,1:binNum4Stat),2),sum(temp_raster(:,beforeEvent:afterEvent),2));
   [all_psth.good(good,:),all_psth.time_vector] = psth(temp_raster,binSize);
   title(['unit_ ' , num2str(good)])
end

if isempty(all_psth.good)
    good = 0;
end


for mua = 1:size(mua_unit_times,1)
   temp_raster = createRaster(event_times,mua_unit_times(mua,:),time_window);
   [all_psth.stats(good+mua,1),all_psth.stats(good+mua,2)] = signrank(sum(temp_raster(:,1:binNum4Stat),2),sum(temp_raster(:,beforeEvent:afterEvent),2));
   all_psth.mua(mua,:) = psth(temp_raster,binSize);
   title(['unit_ ' , num2str(good+mua)])
end

all_psth.United = [all_psth.good;all_psth.mua];



end

