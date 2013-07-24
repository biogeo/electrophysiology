function out = bandlimit(data,sr,band)
% computes band-limited power of a time series
% data is a matrix of data (to be filtered along first non-singleton
% dimension -- i.e., each column a variable)
% sr is the sampling rate (in Hz)
% band is either a named frequency band or a two-element vector of
% frequencies (in Hz)

dt = 1/sr;
fullband = [0.01 120];

if ~exist('band','var')
    band = fullband; %default to full passband
end

if ~ischar(band) % if band is numerical specification
    fband = band; 
else
    switch band
        case 'delta'
            fband = [0.1 4];
        case 'theta'
            fband = [4 8];
        case 'alpha'
            fband = [8 13];
        case 'beta'
            fband = [13 30];
        case 'gamma'
            fband = [30 100];
        otherwise
            fband = fullband;
    end
end
            

%%%%%%%% apply bandpass filter %%%
[b,a]=ellip(2,0.1,40,fband*2/sr);
out = filtfilt(b,a,data); %filter works along first non-singleton dim

