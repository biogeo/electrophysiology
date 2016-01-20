function [stpsd, f, t] = compute_mv_stpsd(x, window, overlap, nfft, Fs, nAvg)
% Calculate a short-time multivariate power spectral density estimate using
% Welch's method.

defaultarg window overlap nfft Fs
defaultarg -value 5 nAvg

[stpdg, f, t] = compute_mv_stperiodogram(x, window, overlap, nfft, Fs);

% First, compensate for using a one-sided periodogram
nf = size(stpdg,2);
if mod(nf,2)
    % Odd number of frequency samples; all points but DC require
    % compensation
    stpdg(:,2:end,:) = 2*stpdg(:,2:end,:);
else
    % Even number of frequency samples; all points but DC and Nyquist
    % require compensation
    stpdg(:,2:end-1,:) = 2*stpdg(:,2:end-1,:);
end

if isempty(Fs)
    stpsd = stpdg / (2*pi);
else
    stpsd = stpdg / Fs;
end

stpsd = sliding_mean(stpsd, nAvg, 3);
t = sliding_mean(t, nAvg); % Dirty and inefficient, can fix later