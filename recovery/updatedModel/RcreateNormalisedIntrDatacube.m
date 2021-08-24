function [amIntrNormcube] = RcreateNormalisedIntrDatacube(amIntrDatacube, normmean, normstd, ninterventions, nmeasures, measures, sigmamethod)

% createNormalisedIntrDatacube - creates the normalised data cube by
% intervention (for each measure)

% for sigma methods 1, 2, & 3 just normalise by mu (as the sigma is
% constant for a given intervention/measure and is incorporated in the
% model objective function
% for sigma methos 4, need to normalise by mu and sigma here as the model
% is using a by day/measure sigma.

amIntrNormcube = amIntrDatacube;
invmeasarray = getInvertedMeasures('BR');

for i = 1:ninterventions
    for m = 1:nmeasures
        if ~ismember(measures.DisplayName(m),invmeasarray)
            if sigmamethod == 4
                amIntrNormcube(i, :, m) = ...
                    (amIntrDatacube(i, :, m) - normmean(i, m)) / normstd(i, m);
            else 
                amIntrNormcube(i, :, m) = ...
                    (amIntrDatacube(i, :, m) - normmean(i, m));
            end
        else
            if sigmamethod == 4
                amIntrNormcube(i, :, m) = ...
                    - (amIntrDatacube(i, :, m) - normmean(i, m)) / normstd(i, m);
            else 
                amIntrNormcube(i, :, m) = ...
                    - (amIntrDatacube(i, :, m) - normmean(i, m));
            end
        end
    end
end

end

