% creates the inputs for the model
% 
% - creates list of interventions with enough data to run model on
% - formats measure data into a datacube of size (npatients, ndays,
% nmeasurestype)
% - outputs this in a .mat file
% 
% Input:
% ------
% clinical and measurements data
% ivandmeasures_recovery_gap*.mat
%
%
% Output:
% -------
% alignmentmodelinputs_recovery_gap*.mat

init;
chosentreatgap = selectTreatmentGap();

[datamatfile, clinicalmatfile, demographicsmatfile] = getRawDataFilenamesForStudy(study);
[physdata, offset] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
[~, cdDrugTherapy, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = loadAndHarmoniseClinVars(clinicalmatfile, subfolder, study);

tic
fprintf('Loading datademographics by patient\n');
load(fullfile(basedir, subfolder, demographicsmatfile));
toc

tic
ivandmeasuresfile = sprintf('%sivandmeasures_recovery_gap%d.mat', study, chosentreatgap);
fprintf('Loading post treatment and measures data\n');
load(fullfile(basedir, subfolder, ivandmeasuresfile));
toc

% useful variables
npatients = max(physdata.SmartCareID);
ndays = max(physdata.ScaledDateNum);

[measures, nmeasures] = createMeasuresTable(physdata);

tic
% create list of interventions with enough data to run model on
fprintf('Creating list of interventions\n');
amInterventions = createListOfInterventionsRecovery(ivandmeasurestable, physdata, offset, data_window);
ninterventions = size(amInterventions,1);
toc

amInterventions = addDrugTherapyInfo(amInterventions, cdDrugTherapy);

tic
% create datacube - 3D array of patients/days/measures for model
% note data cube is 96% empty because the first column contains the ID
% number in range 1:max(ID number) -> no ID = 1 -> very scarce
fprintf('Creating 3D data array\n');
[amDatacube] = createDataCube(physdata, measures, npatients, ndays, nmeasures);
toc

tic
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
outputfilename = sprintf('%salignmentmodelinputs_recovery_gap%d_datawind%d.mat', study, chosentreatgap, data_window );
fprintf('Saving output variables to file %s\n', outputfilename);
save(fullfile(basedir, subfolder,outputfilename), 'amInterventions','amDatacube', 'measures', 'npatients','ndays', 'nmeasures', 'ninterventions');
toc