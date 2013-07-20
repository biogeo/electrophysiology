function P = timefreq(data, freqs, sr, pars)
% given data, a timeseries matrix with each trial a row
% freqs, a vector of frequencies
% sr, the sampling rate of the time series
% pars = [Fb Fc], parameters (see below)
% return a time-frequency matrix representing the average across trials of
% the continuous wavelet transform performed with Morlet wavelets

if ~exist('pars','var')
    pars = [1 1]; % a pretty typical set of parameters
end

if ~exist('freqs','var')
    % for electrophysiology, some common frequency bands:
    % [0.1 4] = delta; [4 8] = theta; [8 13] = alpha; [13 30] = beta;
    % [30 100] = gamma
    freqs = [1:13 15:3:30 35:5:100];
end

% MORLET WAVELET:
% cmor(x) = (pi*Fb)^{-0.5}*exp(2*i*pi*Fc*x)*exp(-(x^2)/Fb) -- from Matlab
Fb = pars(1); %from Cavanagh paper
Fc = pars(2); %from Cavanagh paper
wname = sprintf('cmor%f-%f',Fb,Fc);

scales = sr./freqs; 

% preallocate for speed freqs x times x trials
coefs = nan(length(freqs),size(data,2),size(data,1));
for ind = 1:size(data,1)
    coefs(:,:,ind)=cwt(data(ind,:),scales,wname); %wavelet transform
end

PP = abs(coefs).^2; %calculate power
P = nanmean(PP,3); %average across trials
