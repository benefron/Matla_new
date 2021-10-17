function [all_psth] = mean_psth(event_times,unit_struct,time_window,binSize)

good_unit_times = unit_struct.good.times;
mua_unit_times = unit_struct.mua.times;

all_psth.good = zeros(size(good_unit_times,1),time_window*1000/binSize);
all_psth.mua = zeros(size(mua_unit_times,1),time_window*1000/binSize);
for good = 1:size(good_unit_times,1)
   temp_raster = createRaster(event_times,good_unit_times(good,:),time_window);
   [all_psth.good(good,:),all_psth.time_vector] = psth(temp_raster,10);
end


for mua = 1:size(mua_unit_times,1)
   temp_raster = createRaster(event_times,mua_unit_times(mua,:),time_window);
   all_psth.mua(mua,:) = psth(temp_raster,10);
end




end

