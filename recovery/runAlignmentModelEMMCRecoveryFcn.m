function runAlignmentModelEMMCRecoveryFcn(amRunParameters, interventionslist)

% runs the alignment model (EM version) given a set of run parameters.
%
% This version allows for multiple versions of the latent curves.
% 
% Input:
% ------
% amRunParameters                 parameter file
% *alignmentmodelinputs_gap*.mat  modelinputsmatfile
% *datademographicsbypatient.mat  datademographicsfile
% *dataoutliers.mat               dataoutliersfile
% *LabelledInterventions*.mat     labelledinterventionsfile (manually generated)
% *electivetreatments_gap10.xlsx  electivefile (manually generated)
% 
% Output:
% -------
% mutliple plots
% *obj*.mat                       results (very long file name)

% set the various model run parameters
[mversion, study, treatgap, testlabelmthd, testlabeltxt, ...
    modelinputsmatfile, datademographicsfile, dataoutliersfile, labelledinterventionsfile, electivefile, ...
    sigmamethod, mumethod, curveaveragingmethod, smoothingmethod, datasmoothmethod, ...
    measuresmask, runmode, randomseed, intrmode, modelrun, imputationmode, confidencemode, printpredictions, ...
    offset, align_wind, ~, outprior, heldbackpct, confidencethreshold, nlatentcurves, countthreshold, scenario, vshiftmode, vshiftmax] ...
    = amEMMCSetModelRunParametersFromTableRecovery(amRunParameters);

fprintf('Running Alignment Model %s\n', mversion);
fprintf('\n');

% load the required input data
tic
basedir = setBaseDir();
subfolder = 'MatlabsavedVariables';
fnmodelrun = fullfile(basedir, subfolder, sprintf('%s.mat',modelrun));
fprintf('Loading alignment model Inputs data %s\n', modelinputsmatfile);
load(fullfile(basedir, subfolder, modelinputsmatfile));
fprintf('Loading datademographics by patient %s\n', datademographicsfile);
load(fullfile(basedir, subfolder, datademographicsfile));
fprintf('Loading data outliers %s\n', dataoutliersfile);
load(fullfile(basedir, subfolder, dataoutliersfile));
% if ismember(study, {'SC', 'CL', 'BR'})
%     fprintf('Loading latest labelled test data file %s\n', labelledinterventionsfile);
%     load(fullfile(basedir, subfolder, labelledinterventionsfile), 'amLabelledInterventions');
% end

if ismember(study, {'BR', 'CL'})
    subfolder = sprintf('DataFiles/%s', study);
else
    subfolder = 'DataFiles';
end

fprintf('Loading elective treatment file %s\n', electivefile);
elopts = detectImportOptions(fullfile(basedir, subfolder, electivefile));
elopts.VariableTypes(:, ismember(elopts.VariableNames, {'Hospital'})) = {'char'};
elopts.VariableTypes(:, ismember(elopts.VariableNames, {'PatientNbr', 'ID', 'IVScaledDateNum'})) = {'double'};
amElectiveTreatments = readtable(fullfile(basedir, subfolder, electivefile), elopts);
amElectiveTreatments.ElectiveTreatment(:) = 'Y';
toc

tic
fprintf('Preparing input data\n');

baseplotname = sprintf('%s%s_gp%d_lm%d_sig%d_mu%d_ca%d_sm%d_rm%d_in%d_im%d_cm%d_mm%d_od%d_ou%d_dw%d_nl%d_rs%d_ds%d_ct%d_sc%s_vs%d_vm%.1f', study, mversion, treatgap, testlabelmthd, sigmamethod, mumethod, curveaveragingmethod, ...
    smoothingmethod, runmode, intrmode, imputationmode, confidencemode, measuresmask, offset.down, offset.up, align_wind, nlatentcurves, randomseed, datasmoothmethod, countthreshold, scenario, vshiftmode, vshiftmax);
detaillog = true;

% overwrite ninterventions with interventionslist and amInterventions if
% necessary
if not(isempty(interventionslist))
    amInterventions = amInterventions(interventionslist,:);
    ninterventions = size(amInterventions,1);
end

% 1) select a set of measures, 2) filter corresponding data and 
% 3) compute high level statistics (for each measure)
[amDatacube, measures, nmeasures] = amEMMCPreprocessMeasuresRecovery(amDatacube, amInterventions, measures, ...
    demographicstable, measuresmask, align_wind, npatients, ndays, ninterventions, nmeasures, study);

% create cube for data window data by intervention (for each measure)
[amIntrDatacube] = amEMMCCreateIntrDatacubeRecovery(amDatacube, amInterventions, measures, align_wind, ...
    offset, ninterventions, nmeasures, curveaveragingmethod, datasmoothmethod);

% pre-process intervention table and associated measurement data

% no need to modify max_offset here for Breathe
% 1) (not for BR) calculate data window completeness without interpolated measurements (#actual_points/#max_points*100)
% 2) merge elective treatments into amInterventions
% 3) optionnally remove sequential interventions if intrmode == 2
% Interp measures concerns CLIMB
[amInterventions, amIntrDatacube, ninterventions, intrkeepidx] = amEMMCPreprocessInterventions(amInterventions, ...
    amIntrDatacube, amElectiveTreatments, measures, nmeasures, offset.span, align_wind, ninterventions, intrmode, study);

% populate multiplicative normalisation (sigma) values based on methodology
% selected
normstd = calculateSigmaNormalisation(amInterventions, measures, demographicstable, ninterventions, nmeasures, sigmamethod, study);

% calculate additive normalisation (mu) based on methodology
% and then create normalised data cube.
normmean = calculateMuNormalisationRecovery(amDatacube, amInterventions, measures, demographicstable, ...
    dataoutliers, ninterventions, nmeasures, mumethod, study);

% populate normalised data cube by intervention
[amIntrNormcube] = RcreateNormalisedIntrDatacube(amIntrDatacube, normmean, normstd, ninterventions, nmeasures, measures, sigmamethod);

% populate index array for held back points (to be used for imputation
% held back points is an options triggered by imputationmode = 2, to remove
% heldbackpct (1%) of points prior to finding the latent curve. The
% difference between the value of those points and their corresponding twin
% on the latent curve is an indicator of how well the model worked.
[amHeldBackcube] = RcreateHeldBackcube(amIntrDatacube, offset, align_wind, ninterventions, nmeasures, heldbackpct, imputationmode);

toc
fprintf('\n');

% 1) initialise mean curve arrays (sum, sum squared, count) 
% 2) initialise offset probability distribution
% 3) set the method for multiple latent curve initialisation 
[meancurvesumsq, meancurvesum, meancurvecount, amInterventions, initial_offsets, initial_latentcurve, ...
    animatedmeancurvemean, animatedoffsets, animatedlc, hstg, pdoffset, overall_hist, overall_pdoffset, animated_overall_pdoffset, ...
    vshift, isOutlier, aniterations, run_type] = ...
    RamEMMCInitialiseAlignment(amIntrNormcube, amHeldBackcube, amInterventions, measures, offset, align_wind, ...
    nmeasures, ninterventions, nlatentcurves, runmode, randomseed);

if vshiftmode == 0
    fprintf('Running alignment - without vertical shift\n');
    allowvshift1 = false;
    maxiterations1 = 200;
    maxiterations2 = 0;
elseif vshiftmode == 1
    fprintf('Running alignment - with vertical shift\n');
    allowvshift1 = true;
    maxiterations1 = 200;
    maxiterations2 = 0;
elseif vshiftmode == 2
    fprintf('Running alignment - initially without vertical shift, then with\n');
    allowvshift1 = false;
    maxiterations1 = 200;
    allowvshift2 = true;
    maxiterations2 = 50;
    
end
miniiter = 0;

tic
% align curves, find best offset and best latent
[meancurvesumsq, meancurvesum, meancurvecount, meancurvemean, meancurvestd, amInterventions, ...
    hstg, pdoffset, overall_hist, overall_pdoffset, ...
    animatedmeancurvemean, animatedoffsets, animatedlc, animated_overall_pdoffset, ...
    vshift, isOutlier, pptsstruct, qual, niterations, miniiter] = ...
    amEMMCAlignCurvesRecovery(meancurvesumsq, meancurvesum, meancurvecount, amIntrNormcube, amHeldBackcube, ...
        animatedmeancurvemean, animatedoffsets, animatedlc, animated_overall_pdoffset, ...
        hstg, pdoffset, overall_hist, overall_pdoffset, vshift, isOutlier, ...
        amInterventions, outprior, measures, normstd, offset, align_wind, ...
        nmeasures, ninterventions, nlatentcurves, sigmamethod, smoothingmethod, ...
        runmode, countthreshold, aniterations, maxiterations1, allowvshift1, vshiftmax, miniiter, fnmodelrun);
fprintf('%s - ErrFcn = %.8f\n', run_type, qual);
toc

if maxiterations2 ~= 0
    tic
    [meancurvesumsq, meancurvesum, meancurvecount, meancurvemean, meancurvestd, amInterventions, ...
    hstg, pdoffset, overall_hist, overall_pdoffset, ...
    animatedmeancurvemean, animatedoffsets, animatedlc, animated_overall_pdoffset, ...
    vshift, isOutlier, pptsstruct, qual, niterations2, miniiter] = ...
    amEMMCAlignCurvesRecovery(meancurvesumsq, meancurvesum, meancurvecount, amIntrNormcube, amHeldBackcube, ...
        animatedmeancurvemean, animatedoffsets, animatedlc, animated_overall_pdoffset, ...
        hstg, pdoffset, overall_hist, overall_pdoffset, vshift, isOutlier, ...
        amInterventions, outprior, measures, normstd, offset, align_wind, ...
        nmeasures, ninterventions, nlatentcurves, sigmamethod, smoothingmethod, ...
        runmode, countthreshold, aniterations, maxiterations2, allowvshift2, vshiftmax, miniiter, fnmodelrun);
    niterations = niterations + niterations2;
    fprintf('%s - ErrFcn = %.8f\n', run_type, qual);
    toc
end

% remove points from latent curves that don't have enough underlying data
% points contributing
meancurvemean(meancurvecount  < countthreshold) = nan;
meancurvestd(meancurvecount   < countthreshold) = nan;
meancurvesumsq(meancurvecount < countthreshold) = nan;
meancurvesum(meancurvecount   < countthreshold) = nan;
meancurvecount(meancurvecount < countthreshold) = nan;

% save the zero offset pre-profile to unaligned_profile so all plots show a
% consistent unaligned curve as the pre-profile.
% *** updated to use final curve set assignment with initial uniform offset
% distribution for this ***
% unaligned_profile = profile_pre;
tmp_meancurvesum      = zeros(nlatentcurves, align_wind + offset.span-1, nmeasures);
tmp_meancurvesumsq    = zeros(nlatentcurves, align_wind + offset.span-1, nmeasures);
tmp_meancurvecount    = zeros(nlatentcurves, align_wind + offset.span-1, nmeasures);
tmp_overall_pdoffset  = zeros(nlatentcurves, ninterventions,  offset.span);
for i = 1:ninterventions
    if runmode == 5
        tmp_overall_pdoffset(:, i, :) = 0;
        tmp_overall_pdoffset(amInterventions.LatentCurve(i), i, 1) = 1;
    else
        tmp_overall_pdoffset(amInterventions.LatentCurve(i), i,:) = amEMMCConvertFromLogSpaceAndNormalise(zeros(1, offset.span));
    end
end
for i = 1:ninterventions
    [tmp_meancurvesumsq, tmp_meancurvesum, tmp_meancurvecount] = RamEMMCAddToMean(tmp_meancurvesumsq, tmp_meancurvesum, tmp_meancurvecount, ...
        tmp_overall_pdoffset, amIntrNormcube, amHeldBackcube, zeros(nlatentcurves, ninterventions, nmeasures, offset.span), i, ...
        offset, align_wind, nmeasures, nlatentcurves);
end
[unaligned_profile, ~] = RamEMMCCalcMeanAndStd(tmp_meancurvesumsq, tmp_meancurvesum, tmp_meancurvecount);

plotname = sprintf('%s_obj%.8f', baseplotname, qual);
temp_max_points = zeros(nlatentcurves, 1);
temp_ex_start   = zeros(1, nlatentcurves);

% plot and save aligned curves (pre and post)
amEMMCPlotAndSaveAlignedCurvesRecovery(unaligned_profile, meancurvemean, meancurvecount, meancurvestd, ...
    amInterventions.Offset, amInterventions.LatentCurve, ...
    measures, temp_max_points, offset, align_wind, nmeasures, run_type, temp_ex_start, sigmamethod, plotname, 'Plots', nlatentcurves);

toc
fprintf('\n');

% if ismember(study, {'SC', 'CL', 'BR'})
%     recovery_start = amEMMCCalcExStartsFromTestLabels(amLabelledInterventions(intrkeepidx, :), amInterventions, ...
%                 overall_pdoffset, max_offset, 'Plots', plotname, ninterventions, nlatentcurves);
% else
%    recovery_start = input('Look at best start and enter recovery start: ');
%    fprintf('\n'); 
% end
recovery_start=zeros(nlatentcurves,1);

tic
run_type = 'Best Alignment';
ex_text = sprintf('%d', recovery_start);
plotname = sprintf('%s_ni%d_ex%s_obj%.8f', baseplotname, niterations, ex_text, qual);
plotsubfolder = strcat('Plots', '/', plotname);
mkdir(strcat(basedir, plotsubfolder));
%strcat(measures.ShortName{logical(measures.RawMeas)})

[amInterventions] = RamEMMCCalcConfidenceBounds(overall_pdoffset, amInterventions, offset, ninterventions, confidencethreshold, confidencemode);
[amInterventions] = RamEMMCCalcAbsPredAndBounds(amInterventions, recovery_start, nlatentcurves);
[sorted_interventions, max_points] = amEMMCVisualiseAlignmentDetailRecovery(amIntrNormcube, amHeldBackcube, amInterventions, meancurvemean, ...
    meancurvecount, meancurvestd, overall_pdoffset, measures, offset, align_wind, nmeasures, ninterventions, ...
    run_type, recovery_start, curveaveragingmethod, plotname, plotsubfolder, nlatentcurves);

amEMMCPlotAndSaveAlignedCurvesRecovery(unaligned_profile, meancurvemean, meancurvecount, meancurvestd, ...
    amInterventions.Offset, amInterventions.LatentCurve, ...
    measures, max_points, offset, align_wind, nmeasures, run_type, recovery_start, sigmamethod, plotname, plotsubfolder, nlatentcurves);

% calculate the total number of outliers and the total number of data
% points
[totaloutliers, totalpoints] = RamEMMCCalcTotalOutliers(amIntrDatacube, isOutlier, amHeldBackcube, ...
    amInterventions.Offsetidx, amInterventions.LatentCurve, offset, align_wind, ninterventions);

% calculate imputed probabilities for held back points
[amImputedCube, imputedscore] = RamEMMCCalcImputedProbabilities(amIntrNormcube, amHeldBackcube, ...
    meancurvemean, meancurvestd, normstd, overall_pdoffset, offset, align_wind, ...
    nmeasures, ninterventions, sigmamethod, smoothingmethod, imputationmode, amInterventions.LatentCurve, nlatentcurves);

toc
fprintf('\n');

tic
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
outputfilename = sprintf('%s.mat', plotname);
fprintf('Saving alignment model results to file %s\n', outputfilename);
fprintf('\n');
save(fullfile(basedir, subfolder, outputfilename), 'amDatacube', 'amIntrDatacube', 'amIntrNormcube', ...
    'amHeldBackcube', 'amImputedCube', 'imputedscore', 'amInterventions', 'intrkeepidx', ...
    'meancurvesumsq', 'meancurvesum', 'meancurvecount', 'meancurvemean', 'meancurvestd', 'animatedmeancurvemean', ...
    'initial_offsets', 'initial_latentcurve', 'animatedoffsets', 'animatedlc', 'qual', 'unaligned_profile', ...
    'hstg', 'pdoffset', 'overall_hist', 'overall_pdoffset', 'animated_overall_pdoffset', ...
    'vshift', 'pptsstruct', 'isOutlier', 'outprior', 'totaloutliers', 'totalpoints', ...
    'sorted_interventions', 'max_points', 'normmean', 'normstd', 'measures', 'baseplotname', 'plotname', 'plotsubfolder', ...
    'study', 'mversion', 'treatgap', 'testlabelmthd', 'sigmamethod','offset', 'align_wind', 'recovery_start', 'confidencethreshold', ...
    'sigmamethod', 'mumethod', 'curveaveragingmethod', 'smoothingmethod', 'datasmoothmethod', 'countthreshold', ...
    'measuresmask', 'runmode', 'randomseed', 'intrmode', 'imputationmode', 'heldbackpct', 'confidencemode', 'printpredictions', ...
    'nmeasures', 'ninterventions', 'niterations', 'nlatentcurves', 'scenario', 'vshiftmode', 'vshiftmax', ...
    'modelinputsmatfile', 'datademographicsfile', 'dataoutliersfile', 'labelledinterventionsfile', 'electivefile');
toc
fprintf('\n');

if printpredictions == 1
    tic
    fprintf('Plotting prediction results\n');
    normmode = 1; % plot regular measurement data
    for i=ninterventions
        amEMMCPlotsAndSavePredictions(amInterventions, amIntrDatacube, measures, pdoffset, ...
            overall_pdoffset, hstg, overall_hist, vshift, meancurvemean, normmean, normstd, isOutlier, recovery_start, ...
            i, nmeasures, max_offset, align_wind, sigmamethod, plotname, plotsubfolder, normmode);
    end
    toc
    fprintf('\n');
end

end