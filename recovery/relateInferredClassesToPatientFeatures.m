init;

% load files

[datamatfile, clinicalmatfile, ~] = getRawDataFilenamesForStudy(study);

[cdPatient, cdDrugTherapy, cdMicrobiology, cdAntibiotics, cdAdmissions, ~, cdCRP, ...
    cdClinicVisits, cdOtherVisits, cdEndStudy, cdHghtWght, cdMedications, cdNewMeds, cdUnplannedContact] ...
            = loadAndHarmoniseClinVars(clinicalmatfile, subfolder, study);

load(fullfile(basedir,subfolder,'BRivandmeasures_recovery_gap10.mat'));
load(fullfile(basedir,subfolder,'BRpredictivemodelinputs.mat'), 'pmPatients', 'pmPatientMeasStats');

% load 1 model results file
% specify parameters to spot the right file in dir
nl = 1; dw = 20; mm = 34; vm=0.4; rs=4;
% load all files that are similar except their random seed
dircontent = dir(fullfile(basedir, subfolder, sprintf('*BRvEMMC_gp10_lm1_sig4_mu4_ca2_sm2_rm4_in1_im1_cm2_mm%i_od0_ou15_dw%i_nl%i_rs%i_ds1_ct5_sc22-V_vs1_vm%.1f*.mat',mm,dw, nl, rs, vm)));

%%

if nl == 1
    ModelResultsFile = dircontent(9).name;
    load(fullfile(basedir, subfolder,ModelResultsFile),'amInterventions'); 
    amIntervention_fail = amInterventions;
    amIntervention_fail.LatentCurve(:) = 2;
    
    ModelResultsFile = dircontent(10).name;
    load(fullfile(basedir, subfolder,ModelResultsFile),'amInterventions', 'measures', 'scenario');
    
    amInterventions= [amInterventions; amIntervention_fail];
    amInterventions = sortrows(amInterventions,'IVStartDate','ascend');
    amInterventions = sortrows(amInterventions,'SmartCareID','ascend');
    amInterventions.IntrNbr(:) = 1:size(amInterventions,1);
    
    nl=2;
        
else
    ModelResultsFile = dircontent.name;
    % load file
    load(fullfile(basedir, subfolder,ModelResultsFile),'amInterventions', 'measures', 'scenario'); 
end

% other function inputs
plotname = 'Wilcoxon tests';
plotsubfolder = 'PlotsRecovery/';

%% removing outliers interventions

patientsinit = unique(amInterventions.SmartCareID);
outlier=0;

% remove outlying interventions
if outlier
    outliers = [8 41 55 94 102 104];
    amInterventions(ismember(amInterventions.IntrNbr,outliers),:) = [];
    
    patientsintr = unique(amInterventions.SmartCareID);
    fprintf('Removing %i outlier interventtions, %i out of %i patients remain\n', length(outliers), length(patientsintr), length(patientsinit));
end
patientsintr = unique(amInterventions.SmartCareID);
ninterventions = size(amInterventions,1);


%%

% adding time to treatment response to amInterventions
%rs1
recovery_start_C1 = 6; % defined by expert judgement looking at typical profile
recovery_start_C2 = 2;
% %faillist
recovery_start_C1 = 10; % defined by expert judgement looking at typical profile
recovery_start_C2 = 2;

idxC1 = amInterventions.LatentCurve == 1;
idxC2 = not(idxC1);
amInterventions.TimeToResponse(idxC1) = getTimeToResponse(recovery_start_C1, amInterventions(idxC1,:));
amInterventions.TimeToResponse(idxC2) = getTimeToResponse(recovery_start_C2, amInterventions(idxC2,:));

% from pmPatients, uses only BMI, sex age, equivalent info in
% cdPatient

% from pmpatientsmeasstats, uses only MeasureIndex, {'PatientNbr', 'Study',
% 'ID', 'RobustMax'} and this only for FEV1

% check if all patients form pmpatientsmeasstats are also in
% aminterventions
fprintf(sprintf('%i patients from amInterventions not in pmPatientsMeasStats (should be 0)\n', sum(not(ismember(patientsintr, unique(pmPatientMeasStats.ID)))) ));

wcxpvaltable = amEMMCPlotVariablesVsLatentCurveSetForPaperRecovery(amInterventions, cdPatient, pmPatientMeasStats, ...
    ivandmeasurestable,cdMicrobiology, cdAntibiotics, cdAdmissions, cdCRP, measures, plotname, plotsubfolder,...
    ninterventions, nl, scenario, rs, study); 

function tau = getTimeToResponse(recovery_start, amInterventions)

tau = zeros(size(amInterventions,1),1);
% find idx where offset < recovery start
idxtoreplace = amInterventions.Offset < recovery_start;
tau(idxtoreplace) = recovery_start - amInterventions.Offset(idxtoreplace);

end
