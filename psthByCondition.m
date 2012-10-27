function [psth,t,err] = psthByCondition(spikes,events,condition,minT,maxT,binT)
% psthByCondition  Simple PSTH of neuronal firing by trial condition
% Usage:
%     psth = psthByCondition(spikes,events,condition,minT,maxT,binT)
%     psth = psthByCondition(spikes,events,[],minT,maxT,binT)
%     [psth, t, err] = psthByCondition(...)
% Inputs:
%   spikes: a vector of all spike timestamps for a unit during a session
%   events: a vector of all timestamps for a specific event during a
%       session
%   condition: an m-by-n logical matrix, where m must be equal to
%       numel(events) and n is the number of conditions to evaluate.
%       condition(i,j) iff event(i) is part of condition j.
%   minT: the start of the epoch relative to the event
%   maxT: the end of the epoch relative to the event
%   binT: the bin width for the PSTH.
% Outputs:
%   psth: a T-by-n matrix, where T is the number of bins in the PSTH and n
%       is the number of conditions. psth(:,j) is the PSTH evaluated under
%       condition j.
%   t: a T-by-1 vector of timestamps of the left edges of the bins in psth.
%   err: a T-by-n matrix of the standard error of the mean for each PSTH
%       timepoint and condition.
%
% Example:
%   [psth, t] = psthByCondition(spikes,events,condition,minT,maxT,binT);
%   plot(t,psth);  % Produces a plot of n PSTH lines, one per condition.

if isempty(condition)
    condition = true(numel(events),1);
elseif size(condition,1)~=numel(events)
    error('condition must be m-by-n, where m is number of events')
end
spikes = spikes(:);
events = events(:);
doErr = nargout >= 3;

t = minT:binT:maxT;
psth = zeros(numel(t), size(condition,2));
if doErr
    err = zeros(size(psth));
end
eventPSTH = zeros(numel(t), numel(events));
for k=1:numel(events)
    thisSpikes = spikes(spikes>=minT+events(k) & spikes<=maxT+events(k));
    eventPSTH(:,k) = histc(thisSpikes-events(k), t);
end
eventPSTH = eventPSTH / binT; % Convert to a rate representation
for k=1:size(condition,2)
    psth(:,k) = mean(eventPSTH(:,condition(:,k)'),2);
    if doErr
        err(:,k) = stdErr(eventPSTH(:,condition(:,k)),2);
    end
end
t = t(1:end-1);
psth = psth(1:end-1,:);
if doErr
    err = err(1:end-1,:);
end