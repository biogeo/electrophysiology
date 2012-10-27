function [raster,t] = eventRaster(spikeTS,eventTS,startT,endT,binWidth)

% loopyPSTH   A simple PSTH function using a for loop
% Usage:
%     [psth, binT] = loopyPSTH(spikeTS, eventTS, minT, maxT, binWidth)
%
% Inputs (all values are in seconds):
%     spikeTS:  a vector of spike timestamps
%     eventTS:  a vector of event timestamps
%     minT:     the time, relative to each event, of the start of the PSTH
%               analysis epoch
%     maxT:     the time, relative to each event, of the end of the PSTH
%               analysis epoch
%     binWidth: the width of each bin in the PSTH. Defaults to 1 ms.
% Outputs:
%     raster: a T-by-N matrix of the firing rate, in Hz, of the unit
%           represented in spikeTS, aligned to each timestamp in eventTS,
%           where T is the number of time points in the raster and N is the
%           number of events. mean(raster,2) is a PSTH.
%     t: the time of the left edge of each bin in the raster, in seconds,
%           measured relative to the event represented in eventTS
%
% eventRaster uses the same approach as loopyPSTH, but preserves per-event
% information. Useful for raster plots or for subdividing a session on a
% per-trial basis.
%
% GKA 11/2011

if ~exist('binWidth','var') || isempty(binWidth)
    binWidth = .001;
end

% Get bin left edge times (with one extra bin at maxT):
t = (startT:binWidth:endT)';
% Preallocate a vector for the psth:
raster = zeros(numel(t), numel(eventTS));
for k=1:numel(eventTS)
    % Compute a single-trial firing histogram for the current event
    raster(:,k) = histc(spikeTS, t + eventTS(k));
end
% Convert counts to firing rate:
raster = raster/binWidth;
% Discard extra bin:
t = t(1:end-1);
raster = raster(1:end-1,:);

% As a courtesy to the user, check to make sure binWidth was chosen wisely:
if mod(startT-endT, binWidth)
    warning('eventRaster:fractionalBin', ...
        ['eventRaster: The given binWidth does not fit in the ' ...
        'analysis epoch an integer number of times. The last ' ...
        'incomplete bin is excluded.']);
end