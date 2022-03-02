clc; clear; close all;

tic
studynbr = 4;
study = 'BR';
studyfullname = 'Breathe';
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
[datamatfile, clinicalmatfile, demographicsmatfile] = getRawDataFilenamesForStudy(study);
[physdata, offset] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
[cdPatient, cdDrugTherapy, cdMicrobiology, cdAntibiotics, cdAdmissions, cdPFT, cdCRP, ...
    cdClinicVisits, cdOtherVisits, cdEndStudy, cdHghtWght] = loadAndHarmoniseClinVars(clinicalmatfile, subfolder, study);
toc

tic
% get the date scaling offset for each patient
patientoffsets = getPatientOffsets(physdata);

% extract clinical FEV1 measures and join with offsets to keep only those patients who
% have enough data (ie the patients left after outlier date handling
pclinicalfev = sortrows(cdPFT(:,{'ID', 'LungFunctionDate', 'FEV1'}), {'ID', 'LungFunctionDate'}, 'ascend');
pclinicalfev.Properties.VariableNames{'ID'} = 'SmartCareID';
pclinicalfev = innerjoin(pclinicalfev, patientoffsets);

% create a scaleddatenum to translate the study date to the same normalised
% scale as measurement data scaled date num
pclinicalfev.ScaledDateNum = datenum(pclinicalfev.LungFunctionDate) - offset - pclinicalfev.PatientOffset;

% extract study date and join with offsets to keep only those patients who
% have enough data (ie the patients left after outlier date handling
pstudydate = sortrows(cdPatient(:,{'ID', 'Hospital', 'StudyNumber', 'StudyDate'}), 'ID', 'ascend');
pstudydate.Properties.VariableNames{'ID'} = 'SmartCareID';
pstudydate = innerjoin(patientoffsets, pstudydate);

% create a scaleddatenum to translate the study date to the same normalised
% scale as measurement data scaled date num
pstudydate.ScaledDateNum = datenum(pstudydate.StudyDate) - offset - pstudydate.PatientOffset + 1;

% extract just the weight measures from smartcare data
pmeasuresfev = physdata(ismember(physdata.RecordingType,'FEV1Recording'),{'SmartCareID', 'ScaledDateNum', 'FEV'});
pmeasuresfev.Properties.VariableNames{'FEV'} = 'FEV1';

plotsacross = 3;
plotsdown = 5;
plotsperpage = plotsacross * plotsdown;

subfolder = sprintf('Plots/%s', study);
if ~exist(strcat(basedir, subfolder), 'dir')
    mkdir(strcat(basedir, subfolder));
end
toc

tic
% create plots for all patients
fprintf('FEV Plots for all patients\n');
patientlist = unique(pmeasuresfev.SmartCareID);
filenameprefix = sprintf('%s-CalcClinicalVsHomeFEV1', study);
createAndSaveBreatheFEVPlots(patientlist, pmeasuresfev, pclinicalfev, pstudydate, ...
    plotsacross, plotsdown, plotsperpage, subfolder, filenameprefix);

toc
