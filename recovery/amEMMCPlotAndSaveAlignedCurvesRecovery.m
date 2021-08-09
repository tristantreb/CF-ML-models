function amEMMCPlotAndSaveAlignedCurvesRecovery(profile_pre, meancurvemean, meancurvecount, meancurvestd, offsets, latentcurves, ...
    measures, max_points, offset, align_wind, nmeasures, run_type, ex_start, sigmamethod, plotname, plotsubfolder, nlatentcurves)

% amEMMCPlotAndSaveAlignedCurves - wrapper around the
% amEMPlotAndSaveAlignedCurves to plot for each set of latent curves

for n = 1:nlatentcurves
    tmp_profile_pre    = reshape(profile_pre(n, :, :),    [align_wind + offset.span - 1, nmeasures]);
    tmp_meancurvemean  = reshape(meancurvemean(n, :, :),  [align_wind + offset.span - 1, nmeasures]);
    tmp_meancurvecount = reshape(meancurvecount(n, :, :), [align_wind + offset.span - 1, nmeasures]);
    tmp_meancurvestd   = reshape(meancurvestd(n, :, :),   [align_wind + offset.span - 1, nmeasures]);
    tmp_offsets        = offsets(latentcurves == n);
    
    tmp_plotname = sprintf('%s C%d', plotname, n);
    
    amEMPlotAndSaveAlignedCurvesRecovery(tmp_profile_pre, tmp_meancurvemean, tmp_meancurvecount, tmp_meancurvestd, ...
            tmp_offsets, measures, max_points(n, :), offset, align_wind, nmeasures, run_type, ex_start(n), sigmamethod, tmp_plotname, plotsubfolder);
end

end
