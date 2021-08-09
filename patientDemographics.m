init;

[datamatfile, clinicalmatfile, ~] = getRawDataFilenamesForStudy(study);
[cdPatient, cdDrugTherapy, cdMicrobiology, cdAntibiotics, cdAdmissions, cdPFT, cdCRP, ...
cdClinicVisits, cdOtherVisits, cdEndStudy, cdHghtWght] = loadAndHarmoniseClinVars(clinicalmatfile, subfolder, study);

%% limitation : patients in the clinical data that did not provide any measurement in Breathe

[brphysdata, broffset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);

% list and number of patients with no measurement
N_no_measures= abs(size(cdPatient.ID,1) - size(unique(brphysdata.SmartCareID),1));
patient_no_measures_list = join(string(cdPatient.ID(not(ismember(cdPatient.ID, brphysdata.SmartCareID)))),',');
fprintf('%i patients did not provide any measurement: %s\n', N_no_measures, patient_no_measures_list);
% 
% % update all tables with clean patient list from cdPatient
% cdDrugTherapy = cdDrugTherapy(ismember(cdDrugTherapy.ID, cdPatient.ID),:);
% cdMicrobiology = cdMicrobiology(ismember(cdMicrobiology.ID, cdPatient.ID),:);
% cdAntibiotics = cdAntibiotics(ismember(cdAntibiotics.ID, cdPatient.ID),:);
% cdAdmissions = cdAdmissions(ismember(cdAdmissions.ID, cdPatient.ID),:);
% cdPFT = cdPFT(ismember(cdPFT.ID, cdPatient.ID),:);
% cdCRP = cdCRP(ismember(cdCRP.ID, cdPatient.ID),:);
% cdClinicVisits = cdClinicVisits(ismember(cdClinicVisits.ID, cdPatient.ID),:);
% cdOtherVisits = cdOtherVisits(ismember(cdOtherVisits.ID, cdPatient.ID),:);
% cdEndStudy = cdEndStudy(ismember(cdEndStudy.ID, cdPatient.ID),:);
% cdHghtWght = cdHghtWght(ismember(cdHghtWght.ID, cdPatient.ID),:);

%% 
N = size(unique(cdPatient.ID),1);

% Female
female = sum(ismember(cdPatient.Sex,'Female'));

% Age (yr)
age_mean = mean(cdPatient.CalcAgeExact);
age_std = std(cdPatient.CalcAgeExact);

% BMI (kg/m2
bmi_mean = mean(cdPatient.Weight ./ (cdPatient.Height./100),'omitnan');
bmi_std = std(cdPatient.Weight ./ (cdPatient.Height./100),'omitnan');

% FEV1 (% of predicted)
FEV1PrctPredicted = getMergedFEV1PercentagePredicted(brphysdata, cdPatient, cdPFT);
FEV1_mean = mean(FEV1PrctPredicted.Value,'omitnan');
FEV1_std = std(FEV1PrctPredicted.Value, 'omitnan');
% -> segment based on volume
% <40%
v40andunder = sum(FEV1PrctPredicted.Value<40);
% >=40% to <70%
v4070 = sum(FEV1PrctPredicted.Value>=40 & FEV1PrctPredicted.Value<70);
% >= 70% to <90%
v7090 = sum(FEV1PrctPredicted.Value>=70 & FEV1PrctPredicted.Value<90);
% >= 90%
v90andover = sum(FEV1PrctPredicted.Value>=90);

% Genotype
F508del_homozygous = sum(ismember(cdPatient.CFGene1,'F508del') & ismember(cdPatient.CFGene2,'F508del'));
F508del_heterozygous = sum(ismember(cdPatient.CFGene1,'F508del') | ismember(cdPatient.CFGene2,'F508del')) - F508del_homozygous;
other_mutations =  sum(not(ismember(cdPatient.CFGene1,'F508del')) & not(ismember(cdPatient.CFGene2,'F508del')));
fprintf('\nChecksum on genotype: %i (0 = no error)\n', F508del_homozygous+F508del_heterozygous+other_mutations - N);

% CFQ-R treatment domain score

%% Prescribed CFTR modulators
[Drugsbypatients] = getDrugTherapyInfo(cdDrugTherapy);

symkevi = sum(contains(Drugsbypatients.list,"Symkevi"));
trikafta = sum(contains(Drugsbypatients.list,"Trikafta"));
ivacaftor = sum(contains(Drugsbypatients.list,"Ivacaftor"));
orkambi = sum(contains(Drugsbypatients.list,"Orkambi"));

%% Prescribed medication
 
ab_type = groupsummary(cdAntibiotics,{'ID','Route','HomeIV_s'});

% oral
oral = sum(ismember(ab_type.Route,'Oral'));

% iv
iv = sum(ismember(ab_type.Route,'IV'));
% iv hospital
iv_hosp = sum(ismember(ab_type.Route,'IV') & ismember(ab_type.HomeIV_s,'Yes'));
% home ivs
iv_home = sum(ismember(ab_type.Route,'IV') & ismember(ab_type.HomeIV_s,'No'));
fprintf("\nCheck sum on total IVs: %i (0 = no error)\n", iv_hosp + iv_home - sum(ismember(ab_type.Route,'IV')));

% inhaled
ab_name_route_list = groupsummary(cdAntibiotics,{'ID', 'AntibioticName', 'Route'});
inhaled = sum(ismember(ab_name_route_list.Route,'Oral') ...
    & ( ismember(ab_name_route_list.AntibioticName, 'Colomycin') ...
        | ismember(ab_name_route_list.AntibioticName, 'Colistin') ...
        | ismember(ab_name_route_list.AntibioticName, 'Tobramycin') ));

% azithromycin
ab_list = groupsummary(cdAntibiotics,{'ID', 'AntibioticName'});
azithromycin = sum(ismember(ab_list.AntibioticName,'Azithromycin')); % oral or iv

% list of ab
ab_list_unique = unique(cdAntibiotics.AntibioticName);

% list of injected ab
iv_list_unique = unique(cdAntibiotics.AntibioticName(ismember(cdAntibiotics.Route,'IV')));

% list of orally taken ab
oral_list_unique = unique(cdAntibiotics.AntibioticName(ismember(cdAntibiotics.Route,'Oral')));

%% draw table

alinea = "    ";

r = table('Size',[40 2],...
    'VariableTypes',{'string','string'},...
    'VariableNames',{'Characteristic','Value'});

r = addLine(r,1,"Female",sprintf('%i (%.0f%%)', female, female*100/N));
r = addLine(r,2,"Age (yr)",sprintf('%.1f +/-%.1f', age_mean, age_std));
r = addLine(r,3,"BMI (kg/m2)",sprintf('%.1f +/-%.1f', bmi_mean, bmi_std));

r = addLine(r,4,"FEV1 (% of predicted)",sprintf('%.1f +/-%.1f', FEV1_mean, FEV1_std));
r = addLine(r,5,alinea+"Sub-grouping","");
r = addLine(r,6,alinea+alinea+"< 40%",sprintf('%i (%.0f%%)', v40andunder, v40andunder*100/N));
r = addLine(r,7,alinea+alinea+">= 40% to < 70%",sprintf('%i (%.0f%%)', v4070, v4070*100/N));
r = addLine(r,8,alinea+alinea+">= 70% to < 90%",sprintf('%i (%.0f%%)', v7090, v7090*100/N));
r = addLine(r,9,alinea+alinea+">= 90%",sprintf('%i (%.0f%%)', v90andover, v90andover*100/N));

r = addLine(r,10,"Genotype","");
r = addLine(r,11,alinea+"F508del homozygous",sprintf('%i (%.0f%%)', F508del_homozygous, F508del_homozygous*100/N));
r = addLine(r,12,alinea+"F508del heterozygous",sprintf('%i (%.0f%%)', F508del_heterozygous, F508del_heterozygous*100/N));
r = addLine(r,13,alinea+"Other",sprintf('%i (%.0f%%)', other_mutations, other_mutations*100/N));

r = addLine(r,14,"Prescribed CFTR Modulators","");
r = addLine(r,15,alinea+"Triple Therapy",sprintf('%i (%.0f%%)', trikafta, trikafta*100/N));
r = addLine(r,16,alinea+"Symkevi",sprintf('%i (%.0f%%)', symkevi, symkevi*100/N));
r = addLine(r,17,alinea+"Ivacaftor",sprintf('%i (%.0f%%)', ivacaftor, ivacaftor*100/N));
r = addLine(r,18,alinea+"Okrambi",sprintf('%i (%.0f%%)', orkambi, orkambi*100/N));

r = addLine(r,19,"Prescribed Antibiotics","");
r = addLine(r,20,alinea+"IVs",sprintf('%i (%.0f%%)', iv, iv*100/N));
r = addLine(r,21,alinea+"Inhaled antibiotic",sprintf('%i (%.0f%%)', inhaled, inhaled*100/N));
r = addLine(r,22,alinea+"Azithromycin",sprintf('%i (%.0f%%)', azithromycin, azithromycin*100/N));
% CFQ-R
% dornase alfa
% hypertonic saline
% bronchodilator
% employment

%% save table as excel
filename = sprintf('%sPatientDemographics.xlsx',study);
basedir = setBaseDir();
subfolder = 'ExcelFiles';
writetable(r,fullfile(basedir,subfolder,filename),'Sheet',1,'Range','A1');
fprintf(sprintf('Saved %s as Excel file\n', filename));

%% functions
function FEV1prctpredicted = getMergedFEV1PercentagePredicted(brphysdata, cdPatient, cdPFT)

fprintf('\nFEV1 percentage predicted\n');
% -> calculate PercentagePredicted based on brphysdata
FEV1breathe = getMeasureTable(brphysdata,'FEV1Recording','FEV');
% calculate true value based on mean
func = @(x) mean(x);
FEV1breathe = varfun(func,FEV1breathe,'GroupingVariables','ID','InputVariables','FEV1Recording');
% adds predicted FEV1
FEV1breathe = outerjoin(FEV1breathe, cdPatient, 'Type', 'Left', 'Keys', 'ID', 'RightVariables', {'CalcPredictedFEV1'});
% computes % predicted
FEV1breathe.PercentagePredicted = FEV1breathe.Fun_FEV1Recording ./ FEV1breathe.CalcPredictedFEV1 * 100;

% -> calulate based on clinical data
FEV1clinical = varfun(func,cdPFT,'GroupingVariables','ID','InputVariables','FEV1');
FEV1clinical = outerjoin(FEV1clinical, cdPatient, 'Type', 'Left', 'Keys', 'ID', 'RightVariables', {'CalcPredictedFEV1'});
FEV1clinical.PercentagePredicted = FEV1clinical.Fun_FEV1 ./ FEV1clinical.CalcPredictedFEV1 * 100;

% decide which to use
n_fev1clinical = sum(ismember(cdPatient.ID, FEV1clinical.ID));
n_fev1breathe = sum(ismember(cdPatient.ID, FEV1breathe.ID));
fprintf('Number of patients with FEV1 recording in:\n - the clinical data: %i\n - the breathe data: %i\n', ...
    n_fev1clinical, n_fev1breathe);
n_gain = sum(ismember(cdPatient.ID, FEV1clinical.ID) | ismember(cdPatient.ID, FEV1breathe.ID)) - max(n_fev1clinical,n_fev1breathe);
% currently more values in breathe data, let's use it as a reference and
% add
fprintf('We add %i patients by imputing the missing breathe FEV1 data with clinical FEV1\n', n_gain);
fprintf('Patient lung volume semgentation thus uses %i out of %i\n', n_fev1breathe + n_gain, size(cdPatient,1));
% get indexes from both and values fro clinical
FEV1 = outerjoin(FEV1clinical, FEV1breathe, 'Keys', {'ID'}, 'MergeKeys', 1, 'LeftVariables', {'ID', 'PercentagePredicted'}, 'RightVariables', {'ID', 'PercentagePredicted'});
% impute all breathe not nan indexes into clinical, thereby overwriting
% clinical if applicable
idxtoreplace = isnan(FEV1.PercentagePredicted_FEV1breathe) & ~isnan(FEV1.PercentagePredicted_FEV1clinical);
% 
FEV1.PercentagePredicted_FEV1breathe(idxtoreplace) = FEV1.PercentagePredicted_FEV1clinical(idxtoreplace);
FEV1 = renamevars(FEV1, 'PercentagePredicted_FEV1breathe', 'Value');
FEV1prctpredicted = FEV1(:,{'ID','Value'});

end

function r = addLine(r,i,characteristic, value)
    r.Characteristic(i) = characteristic;
    r.Value(i) = value;
end