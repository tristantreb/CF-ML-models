function [sorted_interventions, max_points] = amEMMCVisualiseAlignmentDetailRecovery(amIntrNormcube, amHeldBackcube, amInterventions, meancurvemean, ...
    meancurvecount, meancurvestd, overall_pdoffset, measures, offset, align_wind, nmeasures, ninterventions, ...
    run_type, ex_start, curveaveragingmethod, plotname, plotsubfolder, nlatentcurves)

% amEMMCVisualiseAlignmentDetail - wrapper around the
% amEMVisualiseAlignmentDetail to plot for each set of latent curves

max_points = zeros(nlatentcurves, align_wind + offset.span-1);
sorted_interventions = struct('Curve', []);

for n = 1:nlatentcurves
    tmp_meancurvemean    = reshape(meancurvemean(n, :, :),    [align_wind + offset.span-1, nmeasures]);
    tmp_meancurvecount   = reshape(meancurvecount(n, :, :),   [align_wind + offset.span-1, nmeasures]);
    tmp_meancurvestd     = reshape(meancurvestd(n, :, :),     [align_wind + offset.span-1, nmeasures]);
    tmp_overall_pdoffset = reshape(overall_pdoffset(n, :, :), [ninterventions, offset.span]);
    tmp_ninterventions   = sum(amInterventions.LatentCurve == n);
    tmp_idx              = amInterventions.LatentCurve == n;
    
    tmp_plotname = sprintf('%s C%d', plotname, n);
    
    if tmp_ninterventions ~= 0 
        [tmp_sorted_interventions, max_points(n, :)] = amEMVisualiseAlignmentDetailRecovery(amIntrNormcube(tmp_idx, :, :), ... 
            amHeldBackcube(tmp_idx, :, :), amInterventions(tmp_idx, :), ...
            tmp_meancurvemean, tmp_meancurvecount, tmp_meancurvestd, tmp_overall_pdoffset(tmp_idx,:), ...
            measures, offset, align_wind, nmeasures, tmp_ninterventions, ...
            run_type, ex_start(n), curveaveragingmethod, tmp_plotname, plotsubfolder);
        
        sorted_interventions(n).Curve = table2array(tmp_sorted_interventions);
    end
    
end

end
