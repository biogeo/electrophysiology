function [stft, f, t] = compute_mv_stft(x, window, overlap, nfft, Fs)
% Compute a multivariate short-time Fourier transform
% x is an array representing a multivariate time-domain signal, where the
% first subscript indexes the univariate components (e.g., channel number),
% and the second indexes time.
% 
% stft is an array representing the short-time Fourier transform (using
% overlapping windows) of each univariate component, such that the first
% subscript indexes the component, the second indexes frequency, and the
% third indexes time.

defaultarg window overlap nfft Fs

nsig = size(x,1);

% Preallocate for univariate signals
uv_stfts = cell(nsig,1);

for i=1:nsig
    [uv_stfts{i}, f, t] = ...
        compute_stft(x(i,:), window, overlap, nfft, Fs);
    uv_stfts{i} = shiftdim(uv_stfts{i}, -1);
end

stft = cat(1, uv_stfts{:});