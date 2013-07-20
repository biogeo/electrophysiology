% This file illustrates a basic use case for the code in the repository. It
% loads lfp and event data from sample files and performs a simple
% time-frequency analysis

%% Load data
% note that all times are already in seconds
load sample_events
load sample_lfp

%% Plot raw data
[T, C] = size(lfp);
sep = 1.5;
lfpmat = bsxfun(@plus, lfp, sep * (1:C));
tstamps = (1:T)/sr;
plot(tstamps,lfpmat)
xlabel('Time (s)')
ylabel('Channel')
set(gca, 'YTick', [])

%% Get data around events
channel = 1;
tseries = lfp(:,channel);

tpre = 1;
tpost = 1.5;

[lfpsplit, binT] = evtsplit(tseries, events, tpre, tpost, sr);

%% Plot
plot(binT, nanmean(lfpsplit), 'k', 'linewidth', 2)
xlabel('Time (s)')
ylabel('Voltage (mV)')

%% Make an averaged spectrogram across trials
freqs = 1:50; % frequencies to analyze
win = 0.5; % length of analysis window (in s)
overlap = 0.95; % percentage window overlap

[S, F, T, P] = avgspecgram(lfpsplit, win, overlap, freqs, sr);

%% Plot
imagesc(binT, F, 10*log10(P))
set(gca, 'ydir', 'normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')

%% Alternately, use a wavelet transform
freqs = [1:13 15:3:30 35:5:50]; % transform takes longer, so use log-spaced frequencies

P = timefreq(lfpsplit, freqs, sr);

%% Plot
pcolor(binT, freqs, 10*log10(P)) %times is a vector of times
shading interp
set(gca,'ydir','normal')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
