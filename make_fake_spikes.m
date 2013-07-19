function data = make_fake_spikes(sr, tottime)
% performs a reality check by feeding white noise through a software
% version of the Plexon MAP filtering; can be used with Offline Sorter,
% automated algorithms, etc. to determine false positive rates

% make data
if ~exist('sr','var')
    sr = 40000; %sampling rate in Hz
end

if ~exist('tottime','var')
    tottime = 10; %sampling time in s
end

sig = randn(1,tottime*sr);

% filter
% make filter
lpf=250; %low cutoff
hpf=4e3; %high cutoff --not what Plexon told me, but what looks like sortclient output

%two ways of doing filter design, equivalent for this case

%low-pass = high-cut
%6-pole Butterworth filter
% [b,a]=butter(6,2*hpf/freq,'low');
% Bf=dfilt.df2(b,a);
[z,p,k]=butter(6,2*hpf/sr,'low'); %normalized frequency is 2f/Fsampling
[sos,g]=zp2sos(z,p,k);
Bf2=dfilt.df2sos(sos,g);

%high-pass = low-cut
%2-pole Linkwitz-Riley = two cascaded 1-pole Butterworth filters
[z,p,k]=butter(1,2*lpf/sr,'high'); %normalized frequency is 2f/Fsampling
[sos,g]=zp2sos(z,p,k);
bb=dfilt.df2sos(sos,g);
LW=dfilt.cascade(bb,bb);

FF = dfilt.cascade(LW,Bf2);

data = filter(FF,sig);

end
