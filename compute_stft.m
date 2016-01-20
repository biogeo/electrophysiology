function [stft, f, t] = compute_stft(x, window, overlap, nfft, Fs)
% Calculate the short-time windowed Fourier transform using spectrogram, as
% a precursor for evaluating short-time Welch's method PSD/CPSD. This is
% essentially a thin wrapper around spectrogram.
%
% N.B., the effective frequency resolution is Fs/length(window), and the
% sampled frequency resolution is Fs/nfft.

% Remove leading singleton dimensions. If x is a row vector, it becomes
% a column vector. If it's higher order than a vector, this will let
% spectrogram do the error handling.

defaultarg window overlap nfft Fs

x = shiftdim(x);

if isempty(window)
    window = hamming(floor(numel(x)/8));
elseif isscalar(window)
    % Window specified only as number of sample points, use a Hamming
    % window of the requested length
    window = hamming(window);
end
window = shiftdim(window);

if isempty(overlap)
    % Default overlap is 50%.
    overlap = 0.5;
end
if overlap < 1
    % Overlap requested as a fraction, convert to number of sample points
    overlap = round(overlap * numel(window));
end

[stft, f, t] = spectrogram(x, window, overlap, nfft, Fs);

% Compensate for window power
winpow = window' * window;
stft = stft ./ sqrt(winpow);

