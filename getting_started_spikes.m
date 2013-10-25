% This file illustrates a basic use case for the code in the repository. It
% loads spike time and event data from sample files and plots a simple
% psth.

%% Load data
% note that all times are already in seconds
load sample_events
load sample_spikes

%% Cut snippets of data before and after the events
tpre = -1; % 1 second before event
tpost = 1.5; % 1.5 seconds after event
dt = 0.001; % 1 ms bins
[psth, binT] = eventRaster(times, events, tpre, tpost, dt);

%% Plot rastergram
% every column of psth is a separate event, every row a time bin
imagesc(binT,1:size(psth,2),psth')
xlabel('Relative Time')
ylabel('Trial')

%% Redo cut with larger bins to lower estimation error
dt = 0.010; % 10 ms bins
[psth, binT] = eventRaster(times, events, tpre, tpost, dt);

%% Calculate average firing rate and standard error
fr = nanmean(psth,2); % mean firing rate across trials
sd = nanstd(psth,[],2); % standard deviation across trials
effsamp = sum(~isnan(psth),2); % number of non-NaN samples at each time point
sem = sd./sqrt(effsamp);

%% Smooth ...
wid = 0.02; % 20 ms smoothing window
frsm = gauss_convolve(fr, wid, dt);
frhi = gauss_convolve(fr + sem, wid, dt);
frlo = gauss_convolve(fr - sem, wid, dt);

%% ... and plot
plot_with_sem(psth', wid, dt, 0, binT);
xlabel('Time (s)')
ylabel('Firing rate (spks/s)')

%% Alternately, use shading
plot_with_sem(psth', wid, dt, 1, binT);
xlabel('Time (s)')
ylabel('Firing rate (spks/s)')