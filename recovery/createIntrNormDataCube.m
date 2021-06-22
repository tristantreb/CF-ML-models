function createIntrNormDataCube(amRunParameters)
% create the array with the normalised data for each intervention
% 
% Input:
% ------
% amRunParameters                 parameter file
% *alignmentmodelinputs_gap*.mat  modelinputsmatfile
% *datademographicsbypatient.mat  datademographicsfile
% *dataoutliers.mat               dataoutliersfile
% *electivetreatments_gap10.xlsx  electivefile (manually generated)
% 
% Output:
% -------
% intrnormdatacube_recovery.mat

% note amDatacube has 10'512'970 values from which 39'089 zeroes,  
% 10'137'509 NaN, 375'461 non NaN

% set the various model run parameters
[mversion, study, treatgap, testlabelmthd, testlabeltxt, ...
    modelinputsmatfile, datademographicsfile, dataoutliersfile, labelledinterventionsfile, electivefile, ...
    sigmamethod, mumethod, curveaveragingmethod, smoothingmethod, datasmoothmethod, ...
    measuresmask, runmode, randomseed, intrmode, modelrun, imputationmode, confidencemode, printpredictions, ...
    max_offset, align_wind, outprior, heldbackpct, confidencethreshold, nlatentcurves, countthreshold, scenario, vshiftmode, vshiftmax] ...
    = amEMMCSetModelRunParametersFromTable(amRunParameters);

fprintf('Running Alignment Model %s\n', mversion);
fprintf('\n');

% load the required input data
tic
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
% fnmodelrun = fullfile(basedir, subfolder, sprintf('%s.mat',modelrun));
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
    

fprintf('Loading elective treatment file %s\n', electivefile);
if ismember(study, {'BR', 'CL'})
    subfolder = sprintf('DataFiles/%s', study);
else
    subfolder = 'DataFiles';
end
elopts = detectImportOptions(fullfile(basedir, subfolder, electivefile));
elopts.VariableTypes(:, ismember(elopts.VariableNames, {'Hospital'})) = {'char'};
elopts.VariableTypes(:, ismember(elopts.VariableNames, {'PatientNbr', 'ID', 'IVScaledDateNum'})) = {'double'};
amElectiveTreatments = readtable(fullfile(basedir, subfolder, electivefile), elopts);
amElectiveTreatments.ElectiveTreatment(:) = 'Y';
toc

tic
fprintf('Preparing input data\n');

baseplotname = sprintf('%s%s_gp%d_lm%d_sig%d_mu%d_ca%d_sm%d_rm%d_in%d_im%d_cm%d_mm%d_mo%d_dw%d_nl%d_rs%d_ds%d_ct%d_sc%s_vs%d_vm%.1f', study, mversion, treatgap, testlabelmthd, sigmamethod, mumethod, curveaveragingmethod, ...
    smoothingmethod, runmode, intrmode, imputationmode, confidencemode, measuresmask, max_offset, align_wind, nlatentcurves, randomseed, datasmoothmethod, countthreshold, scenario, vshiftmode, vshiftmax);
detaillog = true;

% 1) select a set of measures, 2) filter corresponding data and 
% 3) compute high level statistics (for each measure)
[amDatacube, measures, nmeasures] = amEMMCPreprocessMeasuresRecovery(amDatacube, amInterventions, measures, ...
    demographicstable, measuresmask, align_wind, npatients, ndays, ninterventions, nmeasures, study);

% create cube for data window data by intervention (for each measure)
[amIntrDatacube] = amEMMCCreateIntrDatacubeRecovery(amDatacube, amInterventions, measures, align_wind, ...
    max_offset, ninterventions, nmeasures, curveaveragingmethod, datasmoothmethod);

% pre-process intervention table and associated measurement data

% 1) calculate data window completeness without interpolated measurements (#actual_points/#max_points*100)
% 2) optionnally remove sequential interventions if intrmode == 2
% Interp measures concerns CLIMB
[amInterventions, amIntrDatacube, ninterventions, intrkeepidx] = amEMMCPreprocessInterventions(amInterventions, ...
    amIntrDatacube, amElectiveTreatments, measures, nmeasures, max_offset, align_wind, ninterventions, intrmode, study);

% populate multiplicative normalisation (sigma) values based on methodology
% selected
normstd = calculateSigmaNormalisation(amInterventions, measures, demographicstable, ninterventions, nmeasures, sigmamethod, study);

% calculate additive normalisation (mu) based on methodology
% and then create normalised data cube.
% TODO % what happens here?
normmean = calculateMuNormalisationRecovery(amDatacube, amInterventions, measures, demographicstable, ...
    dataoutliers, align_wind, ninterventions, nmeasures, mumethod, study);

% populate normalised data cube by intervention
[amIntrNormcube] = createNormalisedIntrDatacube(amIntrDatacube, normmean, normstd, ...
    max_offset, align_wind, ninterventions, nmeasures, sigmamethod);

tic
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
outputfilename = sprintf('%sintrnormdatacube_recovery.mat', study);
fprintf('Saving Intervention Normalised Data Cube to file %s\n', outputfilename);
fprintf('\n');
save(fullfile(basedir, subfolder, outputfilename), 'amIntrNormcube', 'measures');
toc
end

