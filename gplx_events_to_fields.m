function events = gplx_events_to_fields(events)
% Convenience function for easily converting the "events" field from
% gplx_load to a slighly more usable format, if desired. gplx_load produces
% an events struct with three fields: name, channel, and ts, each having N
% entries, where N is the number of event channels. gplx_events_to_fields
% converts this to a struct with N fields, where each field name is an
% event channel name, and the values are the timestamps (channel number
% information is discarded). This isn't the default behavior of gplx_load
% because the event names may not be guaranteed to be valid Matlab struct
% field names.
% Example
%     data = gplx_load(filename);
%     data.events = gplx_events_to_fields(data.events);

S = [events.name(:)'; events.ts(:)'];
events = struct(S{:});