function [smoothed_db_spl,fft_freq] = compute_corrected_audiogram(filename, sample_freq, window_size, start_time, end_time)
    % Define the reference pressure (20 µPa in Pascals)
    P_ref = 2e-5;

    % Load the microphone response curve
    mic_response = readtable('curve.csv');
    mic_freq = mic_response.Var1 * 1000;  % Convert from kHz to Hz
    mic_db = mic_response.Var2;

    % Load the data
    data = readtable(filename);
    signal = data.Var2;

    % If start_time and end_time are provided, select a subset of the signal
    if nargin >= 4
        start_index = round(start_time * sample_freq) + 1;
        end_index = round(end_time * sample_freq);
        signal = signal(start_index:end_index);
    end

    % Perform Fast Fourier Transform (FFT)
    fft_values = fft(signal);

    % Compute power spectral density (PSD)
    psd = abs(fft_values).^2;

    % Compute the frequencies associated with the PSD values
    fft_freq = (0:length(fft_values)-1) .* sample_freq / length(psd);
    
    % Only keep the positive frequencies (since the spectrum is symmetric)
    half_length = ceil((length(psd)+1) / 2);
    psd = psd(1:half_length);
    fft_freq = fft_freq(1:half_length);

    % Convert PSD to dB SPL
    db_spl = 10 .* log10(psd / (P_ref^2));

    % Correct the dB SPL values for the microphone response
    % Correct the dB SPL values for the microphone response
    mic_response_adjusted = interp1(mic_freq, mic_db, fft_freq, 'linear', 'extrap');
    db_spl_corrected = db_spl - mic_response_adjusted';

    % Smooth the corrected dB SPL values using a moving average
    smoothed_db_spl = movmean(db_spl_corrected, window_size);
    
    % Plot the smoothed audiogram
    % figure;
    % plot(fft_freq, smoothed_db_spl);
    % title(filename);
    % xlabel('Frequency (Hz)');
    % ylabel('dB SPL');
    % xlim([2000, 60000]);
end
