function [amHeldBackcube] = RcreateHeldBackcube(amIntrDatacube, offset, align_wind, ninterventions, nmeasures, heldbackpct, imputationmode)

% createHeldBackcube - creates an index array indicating points to be held
% back for imputation (chosen at random)

% set seed random number generator to a consistent point for
% reproduceability
rng(2);

amHeldBackcube = zeros(ninterventions, align_wind + offset.span-1 + abs(offset.down), nmeasures);

if imputationmode ==2
    for i = 1:ninterventions
        for d = 1:offset.span + align_wind + abs(offset.down)
            for m = 1:nmeasures
                if ~isnan(amIntrDatacube(i, d, m))
                    holdback = rand;
                    if holdback <= heldbackpct
                        amHeldBackcube(i, d, m) = 1;
                    end
                end
            end
        end
    end
end
end

