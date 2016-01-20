function y = sliding_mean(x, n, dim)
% Calculate a simple sliding mean

defaultarg dim;

if isempty(dim)
    dim = find(size(x)>1,1,'first');
end

boxcar = ones(n,1)./n;
y = filter(boxcar, 1, x, [], dim);
y = gad_slicedim(y, dim, n:size(y,dim));

% Turns out this method isn't numerically stable. Apparently convolution
% with a boxcar really is the best method for taking a moving mean in
% Matlab. WTH, Mathworks.

% xx = cumsum(x, dim);
% sz = size(xx);
% sz(dim) = 1;
% pad = zeros(sz);
% xx = cat(dim, pad, xx);
% 
% fullinds = cell(size(sz));
% for i=1:numel(fullinds)
%     fullinds{i} = 1:sz(i);
% end
% indexer = @(ind)[fullinds(1:dim-1), ind, fullinds(dim+1:end)];
% last = indexer(n+1:size(xx,dim));
% first = indexer(1:size(xx,dim)-n);
% y = (xx(last{:}) - xx(first{:})) / n;