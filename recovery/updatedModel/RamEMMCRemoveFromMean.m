function [meancurvesumsq, meancurvesum, meancurvecount] = RamEMMCRemoveFromMean(meancurvesumsq, meancurvesum, meancurvecount, ...
    overall_pdoffset, amIntrCube, amHeldBackcube, vshift, currinter, offset, align_wind, nmeasures, nlatentcurves)

% amEMMCRemoveFromMean - remove an underlying curve from each of the sets of meancurve (sumsq, sum and count) 
% from possible offsets, weighted by the overall probability of each
% offset/meancurve

for n = 1:nlatentcurves
    % place the current intervention curve into every possible offset
    % position, weighted by the probability each offset position is the
    % right one
    
    offsetval = offset.down:offset.up;
    for offsetidx = 1:offset.span
        for i = 1:align_wind + offset.span-1 % - abs(offsetval) <- to be added to avoid pulling in more data from the right
            for m = 1:nmeasures
                 if i - offsetval(offsetidx) > 0 ... % emulates W(k,offset) function
                     && ~isnan(amIntrCube(currinter, i - offsetval(offsetidx), m)) ... % point exists
                     && amHeldBackcube(currinter, i - offsetval(offsetidx), m)==0 % point not held back
                         
                    meancurvesumsq(n, i, m) = meancurvesumsq(n, + i, m) - (((amIntrCube(currinter, i - offsetval(offsetidx), m) + vshift(n, currinter, m, offsetidx)) ^ 2) * overall_pdoffset(n, currinter, offsetidx));
                    meancurvesum(n, i, m)   = meancurvesum(n, i, m)     -  ((amIntrCube(currinter, i - offsetval(offsetidx), m) + vshift(n, currinter, m, offsetidx))      * overall_pdoffset(n, currinter, offsetidx));
                    meancurvecount(n, i, m) = meancurvecount(n, i, m)   - overall_pdoffset(n, currinter, offsetidx); % notesum(overall_pdoffset(n, currinter, :) = 1 if 0 nan values for this measure
                end
            end
        end
    end
end

end
