function result = compute_normalized_audiogram(filename, sample_freq, window_size, segment_length, noise_filename, start_time, end_time)
    % Define the reference pressure (20 µPa in Pascals)
    P_ref = 2e-5;

    % Load the microphone response curve
    mic_response = readtable('curve.csv');
    mic_freq = mic_response.Var1 * 1000;  % Convert from kHz to Hz
    mic_db = mic_response.Var2;

    % Load the ambient noise data
    noise_data = readtable(noise_filename);
    noise_signal = noise_data.Var2;

    % Split the noise signal into segments of the same length as the main signal segments
    segment_samples = round(segment_length * sample_freq);
    num_noise_segments = floor(length(noise_signal) / segment_samples);

    % Compute the dB SPL of the first noise segment to get the length of the dB SPL values
    noise_segment_start = 1;
    noise_segment_end = segment_samples;
    noise_segment = noise_signal(noise_segment_start:noise_segment_end);
    fft_values_noise = fft(noise_segment);
    psd_noise = abs(fft_values_noise).^2;
    half_length_noise = ceil((length(psd_noise)+1) / 2);

    % Initialize the array to store the dB SPL values for each noise segment
    db_spl_noise_segments = zeros(num_noise_segments, half_length_noise);


    for i = 1:num_noise_segments
        % Extract the current noise segment
        noise_segment_start = (i-1)*segment_samples + 1;
        noise_segment_end = i*segment_samples;
        noise_segment = noise_signal(noise_segment_start:noise_segment_end);

        % Compute the dB SPL of the current noise segment
        fft_values_noise = fft(noise_segment);
        psd_noise = abs(fft_values_noise).^2;
        fft_freq_noise = (0:length(psd_noise)-1) .* sample_freq / length(fft_values_noise);
        half_length_noise = ceil((length(psd_noise)+1) / 2);
        psd_noise = psd_noise(1:half_length_noise);
        fft_freq_noise = fft_freq_noise(1:half_length_noise);
        db_spl_noise = 10 .* log10(psd_noise / (P_ref^2));
        mic_response_noise_adjusted = interp1(mic_freq, mic_db, fft_freq_noise, 'linear', 'extrap');
        db_spl_noise_corrected = db_spl_noise - mic_response_noise_adjusted';

        % Store the corrected dB SPL values for this noise segment
        db_spl_noise_segments(i, :) = db_spl_noise_corrected';
    end

    % Compute the average dB SPL of the noise segments
    db_spl_noise_average = mean(db_spl_noise_segments, 1);

    % Load the main data
    data = readtable(filename);
    signal = data.Var2;

    % If start_time and end_time are provided, select a subset of the signal
    if nargin >= 6
        start_index = round(start_time * sample_freq) + 1;
        end_index = round(end_time * sample_freq);
        signal = signal(start_index:end_index);
    end

    % Split the signal into segments
    segment_samples = round(segment_length * sample_freq);
    num_segments = floor(length(signal) / segment_samples);
    db_spl_segments = zeros(num_segments, length(db_spl_noise_corrected));

    for i = 1:num_segments
        % Extract the current segment
        segment_start = (i-1)*segment_samples + 1;
        segment_end = i*segment_samples;
        segment = signal(segment_start:segment_end);

        % Perform Fast Fourier Transform (FFT)
        fft_values = fft(segment, length(fft_values_noise));

        % Compute power spectral density (PSD)
        psd = abs(fft_values).^2;

        % Compute the frequencies associated with the PSD values
        fft_freq = (0:length(psd)-1) .* sample_freq / length(fft_values_noise);
        
        % Only keep the positive frequencies (since the spectrum is symmetric)
        half_length = ceil((length(psd)+1) / 2);
        psd = psd(1:half_length);
        fft_freq = fft_freq(1:half_length);

        % Convert PSD to dB SPL
        db_spl = 10 .* log10(psd / (P_ref^2));

        % Correct the dB SPL values for the microphone response
        mic_response_adjusted = interp1(mic_freq, mic_db, fft_freq, 'linear', 'extrap');
        db_spl_corrected = db_spl - mic_response_adjusted';

        % Subtract the dB SPL of the ambient noise
        db_spl_corrected = db_spl_corrected - db_spl_noise_average';

        % Smooth the corrected dB SPL values using a moving average
        db_spl_smoothed = movmean(db_spl_corrected, window_size);

        % Store the smoothed dB SPL values for this segment
        db_spl_segments(i, :) = db_spl_smoothed;
    end

    % Compute the mean, standard deviation, and 95% confidence interval of the dB SPL values
    mean_db_spl = mean(db_spl_segments, 1);
    std_db_spl = std(db_spl_segments, 0, 1);
    ci_lower = mean_db_spl - 1.96 * std_db_spl / sqrt(num_segments);
    ci_upper = mean_db_spl + 1.96 * std_db_spl / sqrt(num_segments);
    
    % Store the results in a structure
    result.fft_freq = fft_freq;
    result.mean_db_spl = mean_db_spl;
    result.std_db_spl = std_db_spl;
    result.ci_lower = ci_lower;
    result.ci_upper = ci_upper;
    
    % Plot the mean audiogram with confidence intervals
    figure;
    plot(fft_freq, mean_db_spl);
    hold on;
    plot(fft_freq, ci_lower, '--');
    plot(fft_freq, ci_upper, '--');
    title('Mean Audiogram with 95% Confidence Intervals');
    xlabel('Frequency (Hz)');
    ylabel('dB SPL');
    xlim([2000, 60000]);
    legend('Mean Audiogram', 'Lower Confidence Interval', 'Upper Confidence Interval');
end
