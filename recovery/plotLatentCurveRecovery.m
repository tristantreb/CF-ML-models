function [xl, yl] = plotLatentCurveRecovery(ax, max_offset, align_wind, offset, meancurve, xl, yl, colour, linestyle, linewidth, anchor)

% plotLatentCurve - plots the latent curve actual and smoothed

if anchor == 1
    dto = max_offset + align_wind -1;
    dfrom = offset + 1;
else
    dto = max_offset + align_wind - 1 - offset;
    dfrom = 1;
end
    
line(ax, dfrom:dto, ...
    meancurve(1:max_offset + align_wind - 1 - offset), ...
    'Color', colour, ...
    'LineStyle', linestyle, ...
    'LineWidth', linewidth);
                        
%line(ax, ((-1 * (max_offset + align_wind - 1)) + offset: -1), ...
%    (meancurve(1:max_offset + align_wind - 1 - offset)), ...
%    'Color', colour, ...
%    'LineStyle', linestyle, ...
%    'LineWidth', linewidth);

xl = [min(min(-1), xl(1)) max(max((-1 * (max_offset + align_wind - 1))), xl(2))];
xlim(xl);
yl = [min(min(meancurve(1:max_offset + align_wind - 1 - offset) * 0.99), yl(1)) max(max(meancurve(1:max_offset + align_wind - 1 - offset) * 1.01), yl(2))];
if yl(1)==yl(2)
    yl = [yl(1), yl(1)+0.01];
end
ylim(yl);

end

