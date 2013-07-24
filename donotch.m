function out = donotch(data,sr)
% applies notch filtering at 60 and 120 Hz to a time series
% data is a matrix of data (to be filtered along first non-singleton
% dimension -- i.e., each column a variable)
% sr is the sampling rate (in Hz)

dt = 1/sr;

%%%%%%%% apply notch filter %%%%%%%
% notch at 60 Hz
wo = 60/(sr/2); bw = wo/35;
[b,a] = iirnotch(wo,bw);
out = filter(b,a,data)'; %filter works on columns

%notch at 120 Hz
wo = 120/(sr/2); bw = wo/35;
[b,a] = iirnotch(wo,bw);
out = filter(b,a,out)'; %filter works on columns