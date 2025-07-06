function [Y_mag2, half_freqz] = One_sided_fft(Data, fs)
    % Define number of points
    N = length(Data);
    % Define frequency vector
    freqz = linspace(1e-2, fs, N);

    % Perform FFT on data
    Y_mag = fft(Data, N);
    ind = floor(N/2) + 1;
    Y_mag2 = abs(Y_mag(1:ind)) / N;
    Y_mag2(2:ind-1) = 2 * Y_mag2(2:ind-1);
    half_freqz = freqz(1:ind);

    half_freqz = half_freqz(:);
    Y_mag2 = Y_mag2(:);
end