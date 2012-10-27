function [psth,binT] = looplessPSTH(spikeTS,eventTS,minT,maxT,binWidth)
% looplessPSTH   A simple PSTH function without any loops
% Usage:
%     [psth, binT] = looplessPSTH(spikeTS, eventTS, minT, maxT, binWidth)
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
% looplessPSTH uses MATLAB's vectorized operations to compute a PSTH
% without any loops. Although this typically represents an improvement over
% functions with loops, in this case it means replicating the entire
% dataset of spike timestamps once for each event. Despite the high speed
% of MATLAB's vectorized operations, this process means that many
% unnecessary computations are being performed. looplessPSTH is a slight
% improvement over stupidPSTH in terms of speed and readability, but is
% still about 3 times slower than loopyPSTH.
%
% GKA 6/2011

% Reshape spikeTS into a column vector:
spikeTS = spikeTS(:);
% Reshape eventTS into a row vector:
eventTS = eventTS(:)';
% Compute a matrix of spike times relative to each event time:
alignTS = bsxfun(@minus, spikeTS, eventTS);
% Get a vector of the bin left edges, plus one extra bin at maxT:
binT = minT:binWidth:maxT;
% Compute a histogram for our analysis epoch:
psth = histc(alignTS(:)', binT);
% Convert the histogram to a firing rate:
psth = psth/numel(eventTS)/binWidth;
% Discard extra bins:
binT = binT(1:end-1);
psth = psth(1:end-1);
