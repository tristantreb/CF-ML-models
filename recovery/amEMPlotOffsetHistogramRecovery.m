function amEMPlotOffsetHistogramRecovery(ax, offsets)

% amEMPlotOffsetHistogram - plots the histogram of offsets for a given
% alignment model run

ymax = max(groupcounts(offsets)) + 5;
histogram(ax, offsets);
xlim(ax, [min(offsets)-0.5, max(offsets) + 0.5]);
ylim(ax, [0 ymax]);
title(ax, sprintf('Offsets (%d Interventions)', size(offsets, 1)));
set(gca,'FontSize', 6);
end

