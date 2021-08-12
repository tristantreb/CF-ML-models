function [amDatacube, measures, nmeasures] = amEMMCPreprocessMeasuresRecovery(amDatacube, amInterventions, measures, ...
    demographicstable, measuresmask, align_wind, npatients, ndays, ninterventions, nmeasures, study)

% amEMMCPreprocessMeasures - 1) select a set of measures, 2) filter 
% corresponding data and 3) compute high level statistics (for each measure)

% remove temperature readings as insufficient datapoints for a number of
% the interventions
if ismember(study, {'SC', 'TM'})
    idx = ismember(measures.DisplayName, {'Temperature'});
    amDatacube(:,:,measures.Index(idx)) = [];
    measures(idx,:) = [];
    nmeasures = size(measures,1);
    measures.Index = (1:nmeasures)';
end

measures = getMeasuresMask(measuresmask, measures);

% calculate the overall & alignment window std for each measure and store in measures
% table. Also the overall min, max and range values by measure (across all
% patients and days)

for m = 1:nmeasures
    tempdata = [];
    for i = 1:ninterventions
        scid   = amInterventions.SmartCareID(i);
        start = amInterventions.IVScaledDateNum(i);
        periodstart = start;
        tempdata = [tempdata; amDatacube(scid, periodstart:(start + align_wind - 1), m)'];  
    end
    tempdata(1) = [];
    
    measures.AlignWindStd(m) = std(tempdata(~isnan(tempdata)));
    tempdata = reshape(amDatacube(:, :, m), npatients * ndays, 1);
    measures.OverallStd(m) = std(tempdata(~isnan(tempdata)));
    [measures.OverallMin(m), measures.OverallMax(m)] = getMeasureOverallMinMax(demographicstable, measures.Name{m});
    measures.OverallRange(m) = measures.OverallMax(m) - measures.OverallMin(m);
end

end