function amEMPlotOffsetHistogramRecovery(ax, offsets, max_offset)

% amEMPlotOffsetHistogram - plots the histogram of offsets for a given
% alignment model run

nsamples = size(offsets, 1);
ymax = max(groupcounts(offsets)) + 5;
histogram(ax, offsets);
xlim(ax, [-0.5, max_offset - 1 + 0.5]);
ylim(ax, [0 ymax]);
title(ax, sprintf('Offsets (%d Interventions)', size(offsets, 1)));
set(gca,'FontSize', 6);
end

