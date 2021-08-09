function [totaloutliers, totalpoints] = RamEMMCCalcTotalOutliers(amIntrDatacube, isOutlier, amHeldBackcube, offsetsidx, latentcurve, offset, align_wind, ninterventions)

% amEMMCCalcTotalOutliers - sums up the total outliers after alignment
% optimisation has run (as well as the total number of data points) -
% handling multiple sets of latent curves

totaloutliers = 0;
totalpoints   = 0;

for i = 1:ninterventions
       totaloutliers = totaloutliers + sum(sum(isOutlier(latentcurve(i), i, :, :, offsetsidx(i))));
       totalpoints   = totalpoints + sum(sum(~isnan(amIntrDatacube(i, abs(offset.up)+1:abs(offset.up) + align_wind, :))));
end

totalpoints = totalpoints - sum(sum(sum(amHeldBackcube)));

end

