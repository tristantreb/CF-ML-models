function [meancurvesumsq, meancurvesum, meancurvecount, amInterventions, initial_offsets, initial_latentcurve, ...
    animatedmeancurvemean, animatedoffsets, animatedlc, hstg, pdoffset, overall_hist, overall_pdoffset, animated_overall_pdoffset, ...
    vshift, isOutlier, aniterations, run_type] = ...
    RamEMMCInitialiseAlignment(amIntrCube, amHeldBackcube, amInterventions, measures, offset, align_wind, ...
    nmeasures, ninterventions, nlatentcurves, runmode, randomseed)

% amEMMCInitialiseAlignment - function to initialise the various variables
% used in the alignment process

aniterations      = 5000;

meancurvesum      = zeros(nlatentcurves, align_wind + offset.span-1, nmeasures);
meancurvesumsq    = zeros(nlatentcurves, align_wind + offset.span-1, nmeasures);
meancurvecount    = zeros(nlatentcurves, align_wind + offset.span-1, nmeasures);

hstg              = zeros(nlatentcurves, nmeasures, ninterventions, offset.span);
pdoffset          = zeros(nlatentcurves, nmeasures, ninterventions, offset.span);
overall_hist      = zeros(nlatentcurves, ninterventions, offset.span);
overall_pdoffset  = zeros(nlatentcurves, ninterventions, offset.span);

vshift            = zeros(nlatentcurves, ninterventions, nmeasures, offset.span);

animatedmeancurvemean      = zeros(nlatentcurves, align_wind + offset.span-1, nmeasures, aniterations);
animated_overall_pdoffset  = zeros(nlatentcurves, ninterventions, offset.span, aniterations);
animatedoffsets            = zeros(ninterventions, aniterations);
animatedlc                 = zeros(ninterventions, aniterations);

isOutlier         = zeros(nlatentcurves, ninterventions, align_wind, nmeasures, offset.span);

if runmode == 6
    % *** need to initialise amInterventions.LatentCurve and also change how ***
    % *** the pdoffset and overall pdoffset are initialised to include extra ***
    % *** dimension ***
    %run_type = 'Pre-selected Start';
    % save off amInterventions so it doesn't get overwritten
    %tempamintr = amInterventions;
    %load(fnmodelrun, 'overall_hist', 'amInterventions');
    %overall_hstg = overall_hist;
    %for i = 1:ninterventions
    %    overall_pdoffset(i,:) = convertFromLogSpaceAndNormalise(overall_hstg(i,:));
    %end
    %tempamintr.Offset = amInterventions.Offset;
    %tempamintr.LatentCurve = amInterventions.LatentCurve;
    %amInterventions = tempamintr;
elseif runmode == 7 || runmode == 8 || runmode == 10 || runmode == 11
    if runmode == 7
        run_type = 'O-Uniform LC-FEV1Split PM';
        ntiles = nlatentcurves;
    elseif runmode == 8
        run_type = 'O-Uniform LC-Elec_ FEV1Split PM';
        ntiles = nlatentcurves - 1;
    elseif runmode == 10
        run_type = 'O-Uniform LC-FEV1Split SP';
        ntiles = nlatentcurves; 
    elseif runmode == 11
        run_type = 'O-Uniform LC-Elec_ FEV1Split SP';
        ntiles = nlatentcurves - 1;
    end
    fprintf('Creating Upper and Lower 50%% splits for FEV1\n');
    fprintf('Loading Predictive Model Patient Measures Stats\n');
    basedir = setBaseDir();
    subfolder = 'MatlabSavedVariables';
    load(fullfile(basedir, subfolder, 'SCpredictivemodelinputs.mat'), 'pmPatientMeasStats');
    npmpatients = size(unique(pmPatientMeasStats.PatientNbr),1);
    mfev1idx  = measures.Index(ismember(measures.DisplayName, 'LungFunction'));
    fev1max  = pmPatientMeasStats(pmPatientMeasStats.MeasureIndex == mfev1idx, {'PatientNbr', 'Study', 'ID', 'RobustMax'});
    fev1max = sortrows(fev1max, {'RobustMax'}, 'ascend');
    fev1max.NTile(:) = 0;
    for i = 1:npmpatients
        fev1max.NTile(i) = ceil((i * ntiles)/ npmpatients);
    end
    fev1max = sortrows(fev1max, {'Study', 'ID'}, 'ascend');
    lc = innerjoin(amInterventions, fev1max, 'LeftKeys', {'SmartCareID'}, 'RightKeys', {'ID'}, 'LeftVariables', {'SmartCareID'}, 'RightVariables', 'NTile');
    amInterventions.Offset(:) = 0;
    amInterventions.LatentCurve = lc.NTile;
    if runmode == 8 || runmode == 11
        amInterventions.LatentCurve(amInterventions.ElectiveTreatment == 'Y') = ntiles + 1;   
    end    
elseif runmode == 4 || runmode == 5 || runmode == 9   
    % populate pdoffset & overall_pdoffset with uniform prior distribution
    
    % *** need to initialise amInterventions.LatentCurve and also change how ***
    % *** the pdoffset and overall pdoffset are initialised to include extra ***
    % *** dimension ***
    rng(randomseed);
    fprintf('Setting random seed to %d\n', randomseed);
    run_type = 'O-Uniform LC-Random';
    amInterventions.Offset(:) = 0;
    amInterventions.LatentCurve(:) = randi([1, nlatentcurves], [ninterventions, 1]);
else
    fprintf('Unsupported runmode for this version of the alignment model\n');
    return;
end

for i = 1:ninterventions
    if runmode == 9 || runmode == 10 || runmode == 11
        upd = amEMMCConvertFromLogSpaceAndNormalise(zeros(1, offset.span)) / (nlatentcurves + 2);
        for m = 1:nmeasures
            for lc = 1:nlatentcurves
                pdoffset(lc, m, i, :) = upd;
            end
            pdoffset(amInterventions.LatentCurve(i), m, i, :) = upd * 3;   
        end
        for lc = 1:nlatentcurves
            overall_pdoffset(lc, i,:) = upd;
        end
        overall_pdoffset(amInterventions.LatentCurve(i), i, :) = upd * 3;
    else
        for m = 1:nmeasures
            pdoffset(amInterventions.LatentCurve(i), m, i, :) = amEMMCConvertFromLogSpaceAndNormalise(zeros(1, offset.span));
        end
        if runmode == 5
            overall_pdoffset(:, i, :) = 0;
            overall_pdoffset(amInterventions.LatentCurve(i), i, 1) = 1;
        else
            overall_pdoffset(amInterventions.LatentCurve(i), i, :) = amEMMCConvertFromLogSpaceAndNormalise(zeros(1, offset.span));
        end
    end
end

initial_offsets = amInterventions.Offset;
initial_latentcurve = amInterventions.LatentCurve;

animated_overall_pdoffset(:, :, :, 1) = overall_pdoffset;

% calculate initial mean curve over all interventions & prior prob
% distribution for offsets
for i = 1:ninterventions
    [meancurvesumsq, meancurvesum, meancurvecount] = RamEMMCAddToMean(meancurvesumsq, meancurvesum, meancurvecount, ...
        overall_pdoffset, amIntrCube, amHeldBackcube, vshift, i, offset, align_wind, nmeasures, nlatentcurves);
end

end


