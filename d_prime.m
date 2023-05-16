function d = d_prime(hit_rate, false_alarm_rate)
    % Ensure that hit rate and false alarm rate are in the range [0, 1]
    hit_rate = min(max(hit_rate, 0), 1);
    false_alarm_rate = min(max(false_alarm_rate, 0), 1);

    % Calculate the z-score for hit rate and false alarm rate
    z_hit = norminv(hit_rate);
    z_fa = norminv(false_alarm_rate);

    % Calculate d'
    d = z_hit - z_fa;
end
