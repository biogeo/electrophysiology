## electrophysiology
Some minimal Matlab code for doing analysis of electrophysiological datasets.

# Getting Started
Your best bet is probably to take a look at the gists containing typical use cases:
* `getting_started_spikes` shows a simple example of loading spike data and plotting at PSTH.
* `getting_started_LFP` illustrates the process of making a time/frequency plot.

# Loading Plexon Data
For loading data from Plexon files, `gplx_load`, `gplx_events_to_fields`, and `gplx_event_ts_byname` simplify and streamline the process. However, they require that the Plexon [Software Development Kit](http://www.plexon.com/sites/default/files/Plexon%20Offline%20SDKs_0.zip) be installed (and perhaps compiled, for non-Windows users) on your system.

# Sample Data
The files `sample_events`, `sample_spikes`, and `sample_lfp` all contain sample data used to test and run the gists.

# Processing Spikes
Typically, spikes are loaded from files as arrays of timestamps. The functions `eventRaster`, `loopyPSTH` and `looplessPSTH` handle the process of converting these to binned spike counts, usually with the option to select a time window around a given set of events. `psthByCondition` extends this process by returning separate averaged firing rates across task conditions.

# Processing Local Field Potentials
If you need to filter, `donotch` will remove 60Hz (and 120Hz) line noise, and `bandlimit` will filter the signal into a frequency band you specify (as well as accepting arguments like `'alpha'` and `'gamma'` for named frequency bands). 

By analogy with `eventRaster`, `evtsplit` returns a matrix (trials x timebins) of data split around events of interest.

These data can be transformed into plots of power in each frequency band across time by using either the average square of the Fourier Transform (`avgspecgram`) or the similar process for a continuous wavelet transform (`timefreq`). 
