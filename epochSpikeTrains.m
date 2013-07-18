function trains = epochSpikeTrains(spikes, events, window)

trains = cell(size(events));
for i=1:numel(trains)
    trains{i} = spikes(spikes >= events + window(1) & ...
        spikes <= events + window(2));
end