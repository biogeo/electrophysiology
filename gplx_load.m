function plxdata = gplx_load(filename, varargin)
% gplx_load   Load data from a Plexon .plx file
% Usage:
%     data = gplx_load(filename)
%     data = gplx_load(filename, 'Option', value, ...)
% gplx_load is a wrapper for Plexon's Matlab library for reading .plx
% files, which is a suite of functions beginning with plx_*. In many cases,
% it is desirable to load multiple pieces of information from a .plx file,
% which requires multiple, repeated calls to the plx_* functions. This can
% be cumbersome, and prone to introducing bugs, particularly if data are
% being read from multiple files for a single analysis. gplx_load
% simplifies loading data into a single function call, returning a single
% struct including all of the loaded data.
% 
% The default behavior of gplx_load is to load timestamps for all units,
% timestamps for all events, and all A/D continuous channel data, but not
% spike waveforms. This behavior can be customized by supplying any of the
% following options (the default values of which are shown in parentheses):
%     Units ('all'): Which units to load timestamps for. May be any of the
%         following:
%           1. A string or cell array of strings naming the unit(s) to
%              load, where each unit's name is constructed by adding a
%              letter a-z to the channel name, as with Plexon's SortClient
%              and Offline Sorter. E.g., "sig001a" is the first sorted unit
%              on channel sig001.  Append a capital U for the unsorted
%              timestamps ("sig001U").
%           2. An m-by-2 numeric array of [chan, unit] pairs, to load
%              sorted unit number "unit" from spike channel number "chan"
%              (unit 0 is the unsorted timestamps for the channel)
%           3. An m-by-2 cell array of {chan, unit} pairs, to load
%              sorted unit number "unit" from the spike channel named
%              "chan" (unit 0 is again unsorted timestamps). Unit can be a
%              vector in order to load multiple units from one channel.
%           4. The string 'all' to load all units that have at least one
%              recorded timestamp.
%     LoadWaveforms (false): Whether to load spike waveforms as well as
%         timestamps.
%     Events ('all'): Which event channels to load. May be any of the
%         following:
%           1. A string or a cell array of strings naming the event
%              channel(s).
%           2. A numeric vector of channel numbers.
%           3. The string 'all', to load all events. (If the file heppens
%              to include an event 'all', just enclose it in a cell array
%              to load that single file.)
%     ADChannels ('all'): Which continuous data channels to load. Valid
%         are as for 'Events'.
%
% gplx_load produces a struct with the following fields, taken from the
% plx_information function:
%     filename: The name of the loaded file
%     plx_version: The version of the Plexon file
%     wf_freq: The sampling frequency for spike waveforms
%     comment: The text comment entered into the file
%     trodalness: a number indicating single wire, stereotrode, tetrode...
%     points_per_wf: The number of sample points per spike waveform
%     prethreshold_points: The number of sample points prior to threshold
%     spike_peak_v: The peak voltage (mV) of the final spike A/D converter
%     spike_ad_res: Resolution of spike A/D converter, in bits
%     slow_peak_v: The peak voltage (mv) of the analog A/D converter
%     slow_ad_res: Resolution of the analog A/D converter, in bits
%     duration: Duration of the file in seconds
%     datestring: A date/time string for the file
%     datestamp: The file's date/time expressed as a Matlab serial date
% The struct also contains the following fields, which are themselves
% structs describing the file contents:
%     events: Data from event channels:
%         name: 1-by-n cell array of event names
%         channel: 1-by-n numeric array of channel numbers
%         ts: 1-by-n cell array of timestamp vectors for each channel
% 
%     units: Data from spike channels:
%         channel: Channel numbers for each unit
%         channel_unit: The within-channel unit number for each unit
%         names: A cell array of strings fully naming each unit (e.g., 
%                'sig001a')
%         ts: A cell array of spike timestamp vectors for each unit
%       The following fields are present only if LoadWaveforms is true:
%         wave_ts: A cell array of spike timestamp vectors for each unit
%                  identifying times of stored waveforms
%         wave_v: A cell array of N-by-M waveform matrices for N waveforms
%                 sampled at M points, in mV
% 
%     adchans: Data from continuous analog channels, arranged as parallel
%              arrays, each entry holding data for one "chunk" of
%              continuous recording, such that a single channel may appear
%              multiple times (once per chunk):
%         channel: 1-by-n numeric array of channel numbers per chunk
%         name: 1-by-n cell array of channel names per chunk
%         fs: 1-by-n numeric array of A/D sampling frequency per chunk
%         ts: 1-by-n numeric array of chunk start times
%         v: 1-by-n cell array of continuous A/D channel data, in mV.
% 
% Example: Load all spike times, events, and continuous channels
%     data = gplx_load(filename);
%
% Example: Load waveforms as well:
%     data = gplx_load(filename, 'LoadWaveforms', true);
%
% Example: Load timestamps for only unit sig001a
%     data = gplx_load(filename, 'Units', 'sig001a')
%
% Example: Load only a subset of events, and no A/D channels
%     data = gplx_load(filename, 'Events', {'Foo', 'Bar'}, ...
%                      'ADChannels', [])

params = parse_options(varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD HEADER INFORMATION

% Retrieve all "header" fields, the values returned by plx_information
header_fields = {'filename', 'plx_version', 'wf_freq', 'comment', ...
    'trodalness', 'points_per_wf', 'prethreshold_points', ...
    'spike_peak_v', 'spike_ad_res', 'slow_peak_v', 'slow_ad_res', ...
    'duration', 'datestring'};
header_values = cell(size(header_fields));
[header_values{:}] = plx_information(filename);
header = [header_fields; header_values];
% Make room for all "content" fields, so that plxdata always has the same
% fields
content_fields = {'datestamp', 'units', 'events', 'adchan'};
content_values = cell(size(content_fields));
contents = [content_fields; content_values];
% Create plxdata struct, filling in information from the header
plxdata = struct(header{:}, contents{:});
% Set filename to plxdata.filename, in case plx_information chose a new
% file via GUI:
filename = plxdata.filename;
% Convert Plexon date string to more immediately usable Matlab serial date
% number
plxdata.datestamp = datenum(plxdata.datestring, 'mm/dd/yyyy HH:MM:SS');

% Retrieve basic information about DSP (spike), event, and A/D channels
[n_dsp, dsp_chanmap] = plx_chanmap(filename); %#ok<ASGLU>
[n_dsp, dsp_names] = plx_chan_names(filename);
dsp_names = cellstr(dsp_names)';
[n_event, event_chanmap] = plx_event_chanmap(filename); %#ok<ASGLU>
[n_event, event_names] = plx_event_names(filename);
event_names = cellstr(event_names)';
[n_ad, ad_chanmap] = plx_ad_chanmap(filename); %#ok<ASGLU>
[n_ad, ad_names] = plx_adchan_names(filename); %#ok<ASGLU>
ad_names = cellstr(ad_names)';

% Determine number of data points in each channel
[ts_count, wf_count, ev_count, ad_count] = ...
    plx_info(filename, 1); %#ok<ASGLU>
ts_count = ts_count(:,2:end); % Strip off blank column

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD EVENTS DATA

is_event = []; % Will be used to identify missing events from the file
% Get all event indices. After this if block, event_ind should contain the
% index into return values from the plx_* functions, and event_ch should
% contain the actual event channel number to pass to plx_* functions, such
% that all(event_ch == event_chanmap(event_ind)).
if ischar(params.Events)
    if strcmp(params.Events, 'all')
        % Caller requested all events
        event_ind = 1:n_event;
    else
        % Caller requested one event by name
        event_ind = find(strcmp(params.Events, event_names));
        is_event = ~isempty(event_ind);
    end
    event_ch = event_chanmap(event_ind);
elseif iscellstr(params.Events)
    % Caller requested one or more events by name
    [is_event, event_ind] = ismember(params.Events, event_names);
    event_ind = event_ind(is_event);
    event_ch = event_chanmap(event_ind);
else
    % Caller requested one or more events by channel number (presumably; if
    % caller supplied non-numeric input here, errors will occur later)
    event_ch = params.Events(:)';
    [is_event, event_ind] = ismember(event_ch, event_chanmap);
    event_ind = event_ind(is_event);
    event_ch = event_ch(is_event);
end
if ~all(is_event)
    % One or more of the requested events were not actually present in the
    % file. We'll just skip it and raise a warning.
    warning('gplx_load:noSuchEvent', 'Requested event not in file.');
end

% Create events struct
ev_count = ev_count(event_ind);
events.name = event_names(event_ind);
events.channel = event_ch;
events.ts = cell(size(event_ch));
for i=1:numel(event_ch)
    if ev_count(i)
        [n, events.ts{i}] = ...
            plx_event_ts(filename, event_ch(i)); %#ok<ASGLU>
    end
end
plxdata.events = events;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD SPIKE UNIT DATA

% Create unit names matrix from spike channel names and sort letter:
unit_suffixes = cellstr(['U' 'a':'z']');
unit_names = strcat(repmat(dsp_names,27,1),repmat(unit_suffixes,1,n_dsp));
unit_sz = size(unit_names);

is_unit = []; % For each named unit, is name valid?
is_dsp = [];  % For each unit, is DSP channel valid?
% Get all unit indices. After this if block, the following five variables
% should be defined:
%   sort_ind : The 1-based sort index of each unit
%   dsp_ind  : The channel index of each unit in returned Matlab arrays
%   unit_ind : The vector-index for each unit in unit_names
%   sort_num : The 0-based sort number of each unit
%   dsp_num  : The channel number, all(dsp_num == dsp_chanmap(dsp_ind))
if ischar(params.Units)
    if strcmp(params.Units, 'all')
        % Caller requested all non-empty units
        [sort_ind, dsp_ind] = find(ts_count);
    else
        % Caller requested one unit by full name
        [sort_ind, dsp_ind] = find(strcmp(params.Units, unit_names));
    end
    unit_ind = sub2ind(unit_sz, sort_ind, dsp_ind);
    is_unit = ~isempty(unit_ind);
    sort_num = sort_ind - 1;
    dsp_num = dsp_chanmap(dsp_ind);
elseif iscellstr(params.Units)
    % Caller requested one or more units by full name
    [is_unit, unit_ind] = ismember(params.Units, unit_names);
    unit_ind = unit_ind(is_unit);
    [sort_ind, dsp_ind] = ind2sub(unit_sz, unit_ind);
    sort_num = sort_ind - 1;
    dsp_num = dsp_chanmap(dsp_ind);
elseif iscell(params.Units)
    % Caller requested one or more units by {channel name, sort number}
    % pairs
    unit_dsp_names = params.Units(:,1);
    [is_dsp, dsp_ind] = ismember(unit_dsp_names, dsp_names);
    dsp_ind = dsp_ind(is_dsp);
    sort_num_lists = params.Units(:,2);
    n = cellfun(@numel, sort_num_lists);
    dsp_ind = expand_by_counts(dsp_ind, n);
    sort_num = [sort_num_lists{:}];
    sort_ind = sort_num + 1;
    dsp_num = dsp_chanmap(dsp_ind);
    unit_ind = sub2ind(unit_sz, sort_ind, dsp_ind);
else
    % Caller requested one or more units by [channel number, sort number]
    % pairs
    dsp_num = params.Units(:,1);
    sort_num = params.Units(:,2);
    [is_dsp, dsp_ind] = ismember(dsp_num, dsp_chanmap);
    dsp_ind = dsp_ind(is_dsp);
    dsp_num = dsp_num(is_dsp);
    sort_ind = sort_num + 1;
    unit_ind = sub2ind(unit_sz, sort_ind, dsp_ind);
end
is_sort = sort_num <= 26 & sort_num >= 0 & ~mod(sort_num, 1);
sort_num = sort_num(is_sort);
% sort_ind = sort_ind(is_sort); % Not actually used now
% dsp_ind = dsp_ind(is_sort);   % Not actually used now
dsp_num = dsp_num(is_sort);
unit_ind = unit_ind(is_sort);
if ~all(is_dsp)
    % Caller requested a unit from a non-existent DSP channel
    warning('gplx_load:noSuchSpikeChan', ...
        'No such spike channel in file.');
end
if ~all(is_sort)
    % Caller requested a unit by sort number with an invalid sort id
    warning('gplx_load:noSuchSortID', ...
        'Unit sort numbers must be integers between 0 and 26.');
end
if ~all(is_unit)
    % Caller requested a unit by full name that was not present in the file
    warning('gplx_load:noSuchUnit', ...
        'No such unit in file.');
end

% Create units struct
units.channel = dsp_num(:)';
units.channel_unit = sort_num(:)';
units.names = unit_names(unit_ind)';
units.ts = cell(1,numel(dsp_num));
for i=1:numel(units.ts)
    % Retrieve all timestamp data
    [n, units.ts{i}] = ...
        plx_ts(filename, dsp_num(i), sort_num(i)); %#ok<ASGLU>
end
if params.LoadWaveforms
    % Fill in waveform data if requested
    units.wave_ts = cell(size(units.ts));
    units.wave_v = cell(size(units.ts));
    for i=1:numel(units.wave_ts)
        [n, npw, units.wave_ts{i}, units.wave_v{i}] = ...
            plx_waves_v(filename, dsp_num(i), sort_num(i)); %#ok<ASGLU>
    end
end

plxdata.units = units;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD CONTINUOUS A/D CHANNEL DATA

is_ad = [];
% Get all ad channel indices. Following this if block, ad_ind should be the
% index into the ad channel lists, and ad_num should be the actual number
% of the ad channel to supply to plx_* functions.
if ischar(params.ADChannels)
    if strcmp(params.ADChannels, 'all')
        ad_ind = find(ad_count);
    else
        ad_ind = find(strcmp(params.ADChannels, ad_names));
    end
    ad_num = ad_chanmap(ad_ind);
elseif iscellstr(params.ADChannels)
    [is_ad, ad_ind] = ismember(params.ADChannels, ad_names);
    ad_ind = ad_ind(is_ad);
    ad_num = ad_chanmap(ad_ind);
else
    ad_num = params.ADChannels(:)';
    [is_ad, ad_ind] = ismember(ad_num, ad_chanmap);
    ad_ind = ad_ind(is_ad);
    ad_num = ad_num(is_ad);
end
if ~all(is_ad)
    warning('gplx_load:noSuchADChan', ...
        'Requested A/D channel not in file.');
end

ad_fs = zeros(1, numel(ad_num));
ad_ts = cell(1, numel(ad_num));
ad_v  = cell(1, numel(ad_num));
num_chunks = zeros(1, numel(ad_num));
for i=1:numel(ad_num)
    [ad_fs(i), n, ad_ts{i}, ad_fn, all_ad_v] = ...
        plx_ad_v(filename, ad_num(i)); %#ok<ASGLU>
    ad_ts{i} = ad_ts{i}';
    num_chunks(i) = numel(ad_ts{i});
    stop_ind = cumsum(ad_fn);
    start_ind = [1; stop_ind(1:end-1)+1];
    ad_v{i} = cell(1, num_chunks(i));
    for j=1:num_chunks(i)
        ad_v{i}{j} = all_ad_v(start_ind(j):stop_ind(j))';
    end
end
adchan.channel = expand_by_counts(ad_num(:)', num_chunks);
adchan.name = expand_by_counts(ad_names(ad_ind), num_chunks);
adchan.fs = expand_by_counts(ad_fs, num_chunks);
adchan.ts = [ad_ts{:}];
adchan.v  = [ad_v{:}];

plxdata.adchan = adchan;

plx_close(plxdata.filename);


function params = parse_options(args)
defaults = struct( ...
    'LoadWaveforms', false, ...
    'Events', 'all', ...
    'ADChannels', 'all', ...
    'Units', 'all');
% Parse function parameter-value inputs, but try to be as compatible as
% possible for different versions of Matlab/Octave.
if exist('parse_inputs','file')==2
    % Use parse_inputs if available
    params = parse_inputs(args, defaults);
elseif exist('inputParser','file')==2
    % Use inputParser if available
    p = inputParser;
    default_fields = fieldnames(defaults);
    for i=1:numel(default_fields)
        p.addParamValue(default_fields{i}, defaults.(default_fields{i}))
    end
    p.parse(args{:});
    params = p.Results;
else
    % No parameters-by-struct available; just copy p-v pairs.
    params = defaults;
    inputs = reshape(args,2,[]);
    for i=1:size(inputs,2)
        params.(inputs{1,i}) = inputs{2,i};
    end
end


function xx = expand_by_counts(x, n)
% Expand the vector x by a parallel vector n indicating the number of times
% each element should be repeated.
% Examples:
%     expand_by_counts('abcd',[1 2 5 2])
%     ans =
%         abbcccccdd
%     expand_by_counts([5 10 15], [2 1 3])
%     ans =
%         5    5   10   15   15   15
%     expand_by_counts({'foo','bar','baz'},[1 2 3])
%     ans =
%         'foo'    'bar'    'bar'    'baz'    'baz'    'baz'
if isempty(x)
    xx = [];
    return;
end
c = cumsum(n);
d = zeros(1, c(end));
i = [1, c(1:end-1)+1];
d(i) = 1;
xi = cumsum(d);
xx = x(xi);
