function [shortened_seg] = lag_corr(aud_avi,aud_op)

aud_avi_r = resample(double(aud_avi),30000,250000);
time_correct = finddelay(aud_op,aud_avi_r);
time_correct_avi = round(time_correct/30000*250000);
shortened_seg = [aud_avi(time_correct_avi:end);zeros(time_correct_avi-1,1)];
end

