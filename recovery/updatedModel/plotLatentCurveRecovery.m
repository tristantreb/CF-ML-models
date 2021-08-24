function [xl, yl] = plotLatentCurveRecovery(ax, offset, align_wind, meancurve, xl, yl, colour, linestyle, linewidth, anchor)

% plotLatentCurve - plots the latent curve actual and smoothed

if anchor == 1
    dto = align_wind + offset.up;
    dfrom = offset.down + 1;
% else
%     dto = align_wind + offset.span-1;
%     dfrom = 1;
end
    
line(ax, dfrom:dto, ...
    meancurve(1:align_wind + offset.span-1), ...
    'Color', colour, ...
    'LineStyle', linestyle, ...
    'LineWidth', linewidth);
                        
%xl = [min(1, xl(1)) max(max(align_wind + offset.span-1), xl(2))];
xlim(xl);
yl = [min(min(meancurve(1:align_wind + offset.span-1) * 0.99), yl(1)) max(max(meancurve(1:align_wind + offset.span-1) * 1.01), yl(2))];
if yl(1)==yl(2)
    yl = [yl(1), yl(1)+0.01];
end
ylim(yl);

end

