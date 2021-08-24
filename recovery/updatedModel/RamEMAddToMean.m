function [meancurvesumsq, meancurvesum, meancurvecount] = RamEMAddToMean(meancurvesumsq, meancurvesum, meancurvecount, ...
    overall_pdoffset, amIntrCube, amHeldBackcube, currinter, offset, align_wind, nmeasures)

% amEMAddToMean - add an underlying curve to the mean curve (sumsq, sum and count)
% to all possible offsets, weighted by the overall probability of each offset

offsetval = offset.down:offset.up;
for offsetidx = 1:offset.span
    for i = 1:align_wind + offset.span-1 % - abs(offsetval) <- to be added to avoid pulling in more data from the right
        for m = 1:nmeasures
             if i - offsetval(offsetidx) > 0 ... % emulates W(k,offset) function
                 && ~isnan(amIntrCube(currinter, i - offsetval(offsetidx), m)) ... % point exists
                 && amHeldBackcube(currinter, i - offsetval(offsetidx), m)==0 % point not held back

                meancurvesumsq(i, m) = meancurvesumsq(i, m) + (amIntrCube(currinter, i - offsetval(offsetidx), m) ^ 2) * overall_pdoffset(currinter, offsetidx);
                meancurvesum(i, m)   = meancurvesum(i, m)     + amIntrCube(currinter, i - offsetval(offsetidx), m)       * overall_pdoffset(currinter, offsetidx);
                meancurvecount(i, m) = meancurvecount(i, m)   + overall_pdoffset(currinter, offsetidx); % notesum(overall_pdoffset(n, currinter, :) = 1 if 0 nan values for this measure
            end
        end
    end
end

end
