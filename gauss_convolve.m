function y = gauss_convolve(x, sigma, dt)
% Simple gaussian convolution for signal smoothing
% Usage:
%     y = gauss_convolve(x, sigma [, dt])
% Smooth a signal x by convolving it with a Gaussian with standard
% deviation sigma. dt, if supplied, is the sample timestep; if not
% supplied, sigma is assumed to be measured in samples (ie, dt=1).
% gauss_convolve handles the boundaries by replicating the boundary points.

if exist('dt','var')
    % Convert sigma to time steps; otherwise leave it as samples.
    sigma = sigma/dt;
end

% Make the smoothing window 5*sigma in both directions
edge = ceil(5*sigma);
filter = normpdf(-edge:edge,0,sigma);
filter = filter/sum(filter);

% buffer is the amount of edge-padding to perform
buffer = ones(1,edge);

% Get the size of x so we can reshape our result when we're done to respect
% original row or column vectorness
szx = size(x);
x = x(:)';

% Make the padded signal
xx = [buffer*x(1), x, buffer*x(end)];

% Do the convolution, stripping off the buffered points
y = conv(xx, filter, 'valid');
% Restore to original row or column-ness.
y = reshape(y, szx);
