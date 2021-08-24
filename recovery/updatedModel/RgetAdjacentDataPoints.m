function [adjsumsqpt, adjsumpt, adjcountpt, range] = RgetAdjacentDataPoints(meancurvesumsqrow, meancurvesumrow, meancurvecountrow, thispoint, thismeasure, countthreshold, align_wind)

% getAdjacentDataPoints - get and store adjacent data points to a point in the latent
% curve - to be used when there are too few underlying curves contributing

i = 1;
range = 0;
adjsumsqpt = 0;
adjsumpt = 0;
adjcountpt = 0;

% concept:
% idx(i) is odd -> borrow from the right
% idx(i) is even -> borrow from the left

% define boundaries
maxright = align_wind-thispoint;
maxleft = thispoint-1;
idx = 1:2*max(maxright, maxleft);

% create idx sequence: 1 -1 2 -2 3 -3 etc
idx = ceil(idx/2) .* round(double(mod(idx,2) - 0.5));

% remove idx outside align_wind
idx(idx < -maxleft | idx > maxright)=[]; 

while (meancurvecountrow(thispoint) + adjcountpt) < countthreshold

    if i > length(idx)
        i = 1;
    end
    range = idx(i);
    adjsumsqpt = adjsumsqpt + meancurvesumsqrow(thispoint + range);
    adjsumpt   = adjsumpt   + meancurvesumrow(thispoint   + range);
    adjcountpt = adjcountpt + meancurvecountrow(thispoint + range);
    i = i + 1;
end

%fprintf('Retrieved %d adjacent points to point %d for measure %d: %d addtional points\n', range, thispoint, thismeasure, adjcountpt);

end