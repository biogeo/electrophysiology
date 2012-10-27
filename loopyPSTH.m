function [psth, binT] = loopyPSTH(spikeTS, eventTS, minT, maxT, binWidth)
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
%     binWidth: the width of each bin in the PSTH
% Outputs:
%     psth: a vector of the  mean firing rate, in Hz, of the unit
%           represented in spikeTS, aligned to each timestamp in eventTS
%     binT: the time of the left edge of each bin in the PSTH, in seconds,
%           measured relative to the event represented in eventTS
%
% loopyPSTH is an improvement of stupidPSTH. Sliding the epoch window
% rather than the spikes increases the speed by a factor of 5 in a simple
% benchmark with test data. This function is also 3 times faster than
% looplessPSTH.
%
% GKA 6/2011

% Get bin left edge times (with one extra bin at maxT):
binT = minT:binWidth:maxT;
% Preallocate a vector for the psth:
psth = zeros(1,numel(binT));
for k=1:numel(eventTS)
    % Compute a single-trial firing histogram for the current event:
    thisPSTH = histc(spikeTS, binT + eventTS(k));
    % Add it to the running total of counts per bin:
    psth = psth + thisPSTH(:)';
end
% Convert counts to firing rate:
psth = psth/numel(eventTS)/binWidth;
% Discard extra bin:
binT = binT(1:end-1);
psth = psth(1:end-1);

% As a courtesy to the user, check to make sure binWidth was chosen wisely:
if mod(maxT-minT, binWidth)
    warning('loopyPSTH:fractionalBin', ...
        ['loopyPSTH: The given binWidth does not fit in the ' ...
        'analysis epoch an integer number of times. The last ' ...
        'incomplete bin is excluded.']);
end