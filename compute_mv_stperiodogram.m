function [stpdg, f, t] = compute_mv_stperiodogram(x, window, overlap, nfft, Fs)
% Compute a full short-time periodogram from a multivariate signal,
% including cross-periodograms.
% 
% x is an array representing the multivariate signal, with the first
% subscript indexing univariate components, and the second indexing time.
% 
% stpdg is an array representing the full periodogram, including
% cross-periodograms, of the multivariate signal. For memory efficiency,
% only the unidirectional pair-wise comparisons are calculated. That is,
% the cross-periodogram of <x(2,:), x(1,:)> is computed, but not <x(1,:),
% x(2,:)>, the latter being merely the complex conjugate of the former.
%
% The first subscript of stpdg indexes the pairwise comparison, such that
% if x is N-by-T, and [I,J] = find(tril(ones(N))), then stpdg(k,:,:) is the
% cross-periodogram of <x(I(k),:), x(J(k),:)>. 

defaultarg window overlap nfft Fs;

[stft, f, t] = compute_mv_stft(x, window, overlap, nfft, Fs);

sz = size(stft);
nsig = sz(1);
npairs = nsig * (nsig+1) / 2;
sz(1) = npairs;

stpdg = zeros(sz);

[I, J] = find(tril(ones(nsig)));

for k=1:npairs
    sig1 = stft(I(k),:,:);
    sig2 = stft(J(k),:,:);
    stpdg(k,:,:) = sig1 .* conj(sig2);
end
