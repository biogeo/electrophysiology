function [S,F,T,P] = avgspecgram(data,window,percover,freqs,sr)
% computes averaged spectrogram across rows of data
% data is a trials x time matrix
% window is a time window over which to do spectrogram (in s)
% percover is the percent overlap between successive windows
% freqs is a vector of frequencies at which to evaluate the FFT
% sr is the sampling rate (in Hz)
% outputs: S, F, T, and P are defined by the outputs from matlab's native
% spectrogram function, save that they are averages (in the case of S and P)

% calculate some handy numbers:
winlen = ceil(window * sr);
noverlap = floor(winlen * percover);

% do something a little dirty and just throw away incomplete rows
hasnans = logical(sum(isnan(data),2));
data(hasnans,:) = [];

% now loop over rows
for ind = 1:size(data,1)
   dd = data(ind,:);
   [ss{ind},F,T,pp{ind}] = spectrogram(dd,winlen,noverlap,freqs,sr);
end

smat = cat(3,ss{:}); %stack all matrices on top of each other
pmat = cat(3,pp{:});
S = mean(smat,3);
P = mean(pmat,3);

