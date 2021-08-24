function [meancurvemean, meancurvestd] = RamEMMCCalcMeanAndStd(meancurvesumsq, meancurvesum, meancurvecount)

% amEMMCCalcMeanAndStd - recalc meancurvemean and meancurvestd arrays

% recalculate mean curve and std by day
meancurvemean = meancurvesum ./ meancurvecount;
meancurvestd  = (abs((meancurvesumsq ./ meancurvecount) - (meancurvemean .* meancurvemean))) .^ 0.5;

end
