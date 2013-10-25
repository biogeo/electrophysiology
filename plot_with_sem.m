function plot_with_sem(x, smwid, dt, type, binT)
% plots the data in x with lines or a shaded patch representing +/- 1
% standard error of the mean
% x is a trials x timepoints matrix of data; the function will plot 
% mean(x) +/- std(x)/sqrt(N)
% smwid is a smoothing width in absolute units (smwid = 0 => no smoothing)
% dt is the time step between successive samples
% type = 0 for dashed lines, 1 for shaded patch
% binT is the time axis (optional)

[ntrials, npts] = size(x);

xm = nanmean(x, 1);
sd = nanstd(x, [], 1);
effsamp = sum(~isnan(x), 1);
sem = sd./sqrt(effsamp);

% smooth
if smwid ~= 0
    xsm = gauss_convolve(xm, smwid, dt);
    xhi = gauss_convolve(xm + sem, smwid, dt);
    xlo = gauss_convolve(xm - sem, smwid, dt);
else
    xsm = xm;
    xhi = xm + sem;
    xlo = xm - sem;
end

if exist('binT', 'var')
    binT = binT(:)';  % make sure binT is a row vector
else
    binT = 1:size(xm, 2);
end

% now plot
hold all
if type == 1
    Xptch = [binT fliplr(binT)];
    Yptch = [xlo fliplr(xhi)];

    if ntrials > 1
        patch(Xptch, Yptch, [0 0 0], 'Facealpha', ...
            0.25, 'EdgeColor', 'none'); %black patch with 25% opacity and no border
    end
    plot(binT, xsm, 'k', 'linewidth',2)
    
else
    plot(binT, xsm, 'k', 'linewidth', 2);
    plot(binT, xhi, 'color', 0.5*[1 1 1], 'linestyle', '--', 'linewidth',1)
    plot(binT, xlo, 'color', 0.5*[1 1 1], 'linestyle', '--', 'linewidth',1)
end
hold off