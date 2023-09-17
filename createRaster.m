function [Raster_mat,event_matrix] = createRaster(event_times,unit_times,time_window,eventPlot)

time_vector = [1:time_window*30000] - (time_window*15000);
event_matrix = ones(length(event_times),length(time_vector));
event_matrix = event_matrix.*time_vector;
try
    event_matrix = event_matrix + event_times;
catch
    event_matrix = event_matrix + event_times';
end

Raster_mat = ismember(event_matrix,unit_times);


y_multiple = 1:length(event_times);
for_plot = Raster_mat.*y_multiple';
figure
hold on
time_vector = time_vector/30000;
for event=1:eventPlot
    scatter(time_vector(for_plot(event,:)>0),nonzeros(for_plot(event,:)),'.','k')

end

end

