function ts = gplx_event_ts_byname(filename, events)
% Retrieve event timestamps from a .plx file using the events' names.
% Usage:
%     ts = gplx_event_ts_byname(filename, events)
% This is a convenience function to make retrieving event timestamps by
% name instead of by event channel number slightly easier.
%     filename is the name of the .plx file to read
%     events can be either a string or a cell array of strings, specifying
%            the event(s) to read
%     ts is either an n-by-1 array of timestamps (if events is a string) or
%            an m-by-1 cell array of such timestamp arrays (if events is a
%            cell array of strings).

[n, event_names]   = plx_event_names(filename);
event_names = cellstr(event_names);
[n, event_chanmap] = plx_event_chanmap(filename);
if iscellstr(events)
    [in_file, event_id] = ismember(events, event_names);
    ts = cell(numel(events),1);
    for i=1:numel(events)
        if ~in_file(i)
            warning('gplx:NoSuchEvent', 'No such event in PLX file');
            continue;
        end
        [n, ts{i}] = plx_event_ts(filename, event_chanmap(event_id(i)));
    end
else
    event_id = find(strcmp(events, event_names));
    if isempty(event_id)
        warning('gplx:NoSuchEvent', 'No such event in PLX file');
        ts = [];
    else
        [n, ts] = plx_event_ts(filename, event_chanmap(event_id));
    end
end
