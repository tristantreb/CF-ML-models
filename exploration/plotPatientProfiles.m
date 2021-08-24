% plot a profile of each measures for each patient over the whole enrolment period
% 
% script extracted and adapted from visualisePatientSummary.m
% 
% Input:
% ------
% clinical and home measurements data
% 
% Output: one plot per patient (containing all measures)
% -------

init;

[studynbr, study, ~] = selectStudy();
[datamatfile, clinicalmatfile, demographicsmatfile] = getRawDataFilenamesForStudy(study);
[physdata, offset, physdata_predateoutlierhandling] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
[cdPatient, cdDrugTherapy, cdMicrobiology, cdAntibiotics, cdAdmissions, cdPFT, cdCRP, ...
    cdClinicVisits, cdOtherVisits, cdEndStudy, cdHghtWght] = loadAndHarmoniseClinVars(clinicalmatfile, subfolder, study);

% Tristan's function to harmonise drug therapy namings - temporary until 
% REDcap is active
cdDrugTherapy.DrugTherapyType = cleanDrugTherapyNamings(cdDrugTherapy.DrugTherapyType);

tic
fprintf('Loading demographic data by patient\n');
load(fullfile(basedir, subfolder, demographicsmatfile), 'demographicstable', 'overalltable');
toc
fprintf('\n');

subfolder = sprintf('Plots/%s', study);
if ~exist(strcat(basedir, subfolder), 'dir')
    mkdir(strcat(basedir, subfolder));
end

runmode = input('Which patients to run for 1) Those with enough data 2) Those without enough data ?');
if runmode ~= 1 & runmode ~= 2
    fprintf('Invalid entry')
    return;
end
if runmode == 1
    patientlist = unique(physdata.SmartCareID);
elseif runmode == 2
    goodpatients = unique(physdata.SmartCareID);
    patientlist = unique(physdata_predateoutlierhandling.SmartCareID(~ismember(physdata_predateoutlierhandling.SmartCareID, goodpatients)));
    physdata = physdata_predateoutlierhandling;
end

patientoffsets = getPatientOffsets(physdata);
cvcol   = [0.94  0.52  0.15];
admcol  = [0.694 0.627 0.78]; 
ivcol   = [1     0     0   ];
oralcol = [1     0.85  0   ];
trplcol = [0     1     0   ];
drugcol = [0     0.8   0.6 ];

for i = 1:size(patientlist,1)

    tic
    scid       = patientlist(i);
    fprintf('Creating patient summary for patient %d\n', scid);
     poffset    = patientoffsets.PatientOffset(patientoffsets.SmartCareID == scid);
     studyid   = cdPatient.StudyNumber{cdPatient.ID == scid};
     hospital   = cdPatient.Hospital{cdPatient.ID == scid};
%     sex        = cdPatient.Sex{cdPatient.ID == scid};
     spstart    = cdPatient.StudyDate(cdPatient.ID == scid);
     spstartdn  = datenum(spstart) - offset - poffset + 1;
%     spend      = cdPatient.StudyDate(cdPatient.ID == scid)+days(183);
     spenddn    = spstartdn + 183;
%     if studynbr == 3 || studynbr == 4
%         spendstatus = ' ';
%     else
%         spendstatus = cdEndStudy.EndOfStudyReason{cdEndStudy.ID == scid};
%     end
%     hmstart    = min(physdata.Date_TimeRecorded(physdata.SmartCareID == scid));
     hmstartdn  = min(physdata.ScaledDateNum(physdata.SmartCareID == scid));
%     hmend      = max(physdata.Date_TimeRecorded(physdata.SmartCareID == scid));
     hmenddn    = max(physdata.ScaledDateNum(physdata.SmartCareID == scid));
%     
%     spclinicvisits  = size(cdClinicVisits.ID(cdClinicVisits.ID == scid & cdClinicVisits.AttendanceDate >= spstart & cdClinicVisits.AttendanceDate <= spend),1);
%     sphospadm       = size(cdAdmissions.ID(cdAdmissions.ID == scid & cdAdmissions.Admitted >= spstart & cdAdmissions.Admitted <= spend),1);
%     spivab          = size(cdAntibiotics.ID(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, 'IV') & cdAntibiotics.StartDate >= spstart & cdAntibiotics.StartDate <= spend),1);
%     sporalab        = size(cdAntibiotics.ID(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, 'Oral') & cdAntibiotics.StartDate >= spstart & cdAntibiotics.StartDate <= spend),1);
%     spcpftmeas      = size(cdPFT.ID(cdPFT.ID == scid & cdPFT.LungFunctionDate >= spstart & cdPFT.LungFunctionDate <= spend),1);
%     spccrpmeas      = size(cdCRP.ID(cdCRP.ID == scid & cdCRP.CRPDate >= spstart & cdCRP.CRPDate <= spend),1);
%     
%     hmclinicvisits  = size(cdClinicVisits.ID(cdClinicVisits.ID == scid & cdClinicVisits.AttendanceDate >= hmstart & cdClinicVisits.AttendanceDate <= hmend),1);
%     hmhospadm       = size(cdAdmissions.ID(cdAdmissions.ID == scid & cdAdmissions.Admitted >= hmstart & cdAdmissions.Admitted <= hmend),1);
%     hmivab          = size(cdAntibiotics.ID(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, 'IV') & cdAntibiotics.StartDate >= hmstart & cdAntibiotics.StartDate <= hmend),1);
%     hmoralab        = size(cdAntibiotics.ID(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, 'Oral') & cdAntibiotics.StartDate >= hmstart & cdAntibiotics.StartDate <= hmend),1);
%     hmcpftmeas      = size(cdPFT.ID(cdPFT.ID == scid & cdPFT.LungFunctionDate >= hmstart & cdPFT.LungFunctionDate <= hmend),1);
%     hmccrpmeas      = size(cdCRP.ID(cdCRP.ID == scid & cdCRP.CRPDate >= hmstart & cdCRP.CRPDate <= hmend),1);
%     
%     allclinicvisits = size(cdClinicVisits.ID(cdClinicVisits.ID == scid),1);
%     allhospadm      = size(cdAdmissions.ID(cdAdmissions.ID == scid),1);
%     allcpftmeas     = size(cdPFT.ID(cdPFT.ID == scid),1);
%     allccrpmeas     = size(cdCRP.ID(cdCRP.ID == scid),1);
%     allivab         = size(cdAntibiotics.ID(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, 'IV')),1);
%     alloralab       = size(cdAntibiotics.ID(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, 'Oral')),1);
%     
%     allcavgfev1     = mean(cdPFT.FEV1(cdPFT.ID == scid));
%     allcstdfev1     = std(cdPFT.FEV1(cdPFT.ID == scid));
%     allcminfev1     = min(cdPFT.FEV1(cdPFT.ID == scid));
%     allcmaxfev1     = max(cdPFT.FEV1(cdPFT.ID == scid));
%     
%     allcavgfev1_    = mean(cdPFT.CalcFEV1_(cdPFT.ID == scid));
%     allcstdfev1_    = std (cdPFT.CalcFEV1_(cdPFT.ID == scid));
%     allcminfev1_    = min(cdPFT.CalcFEV1_(cdPFT.ID == scid));
%     allcmaxfev1_    = max(cdPFT.CalcFEV1_(cdPFT.ID == scid));
%     
%     allcavgcrp      = mean(cdCRP.NumericLevel(cdCRP.ID == scid));
%     allcstdcrp      = std(cdCRP.NumericLevel(cdCRP.ID == scid));
%     allcmincrp      = min(cdCRP.NumericLevel(cdCRP.ID == scid));
%     allcmaxcrp      = max(cdCRP.NumericLevel(cdCRP.ID == scid));
%     
%     hmtotal         = size(physdata.SmartCareID(physdata.SmartCareID == scid),1);
%     hmduration      = max(physdata.DateNum(physdata.SmartCareID == scid)) - min(physdata.DateNum(physdata.SmartCareID == scid));
%     hmavgperday     = hmtotal/hmduration;
                                                
    microbiology = unique(cdMicrobiology.Microbiology(cdMicrobiology.ID==scid));
           
    mplotsacross = 1;
    mplotsdown = 9;
    mplotsperpage = mplotsacross * mplotsdown;
    
    measures = [{'HasColdOrFluRecording'};{'CoughRecording'};{'WellnessRecording'};{'FEV1Recording'};{'O2SaturationRecording'};{'RestingHRRecording'};{'WeightRecording'};{'TemperatureRecording'};{'CalorieRecording'}];
    npages = ceil(size(measures, 1) / mplotsperpage);
    page = 2;
    filenameprefix = sprintf('%s-Patient Summary - ID %d (%s) Hosp %s', study, scid, studyid, hospital);
    imagefilename = sprintf(filenameprefix);
    
    daysfrom = min(spstartdn, hmstartdn) - 14;
    daysto   = max(spenddn, hmenddn) + 14;
    xl = [daysfrom daysto];
    
    cvset     = cdClinicVisits(cdClinicVisits.ID == scid,:);
    cvset.AttendanceDatedn = datenum(cvset.AttendanceDate) - offset - poffset + 1;
    admset    = cdAdmissions(cdAdmissions.ID == scid,:);
    admset.Admitteddn = datenum(admset.Admitted) - offset - poffset + 1;
    admset.Dischargedn = datenum(admset.Discharge) - offset - poffset + 1;
    admdates = unique([admset.Admitteddn ; admset.Dischargedn]);
    ivabset   = cdAntibiotics(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, {'IV'}),:);
    ivabset.Startdn = datenum(ivabset.StartDate) - offset - poffset + 1;
    ivabset.Stopdn = datenum(ivabset.StopDate) - offset - poffset + 1;
    oralabset = cdAntibiotics(cdAntibiotics.ID == scid & ismember(cdAntibiotics.Route, {'Oral'}),:);
    oralabset.Startdn = datenum(oralabset.StartDate) - offset - poffset + 1;
    oralabset.Stopdn = datenum(oralabset.StopDate) - offset - poffset + 1;
    trplset = cdDrugTherapy(cdDrugTherapy.ID == scid & ismember(cdDrugTherapy.DrugTherapyType, {'Triple Therapy'}),:);
    trplset.Startdn = datenum(trplset.DrugTherapyStartDate) - offset - poffset + 1;
    drugset = cdDrugTherapy(cdDrugTherapy.ID == scid & ~ismember(cdDrugTherapy.DrugTherapyType, {'Triple Therapy'}),:);
    drugset.Startdn = datenum(drugset.DrugTherapyStartDate) - offset - poffset + 1;
    
    f = figure('Name',filenameprefix,'DefaultAxesFontSize',12,'Position', [1 1 2000 1000]);

    % get all measures so the plots for each appear in a consistent place
    % across all patients
    for m = 1:size(measures,1)
        measure = measures{m};
        column = getColumnForMeasure(measure);
        scdata = physdata(physdata.SmartCareID == scid & ismember(physdata.RecordingType, measure), :);
        scdata = scdata(:, {'SmartCareID','ScaledDateNum' 'Date_TimeRecorded', column});
        scdata.Properties.VariableNames{column} = 'Measurement';
        
        if size(scdata,1) > 0
            ax = subplot(mplotsdown, mplotsacross, m - (page-2) * mplotsperpage,'Parent',f);
            hold on;
            xlim(xl);
            rangelimit = setMinYDisplayRangeForMeasure(measure);
            yl = setYDisplayRange(min(scdata.Measurement), max(scdata.Measurement), rangelimit);
            ylim(yl);
            ylabel(ax, replace(measure, 'Recording', ''));
            
            plot(ax, scdata.ScaledDateNum, scdata.Measurement, ...
                'Color', [0, 0.65, 1], ...
                'LineStyle', ':', ...
                'Marker', 'o', ...
                'LineWidth',1,...
                'MarkerSize',2,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','g');
    
            for a = 1:size(ivabset,1)
                fill(ax, [ivabset.Startdn(a) ivabset.Stopdn(a) ivabset.Stopdn(a) ivabset.Startdn(a)], ...
                          [yl(1) yl(1) yl(2) yl(2)], ivcol, 'FaceAlpha', '0.1', 'EdgeColor', 'none');
            end
            for a = 1:size(oralabset,1)
                fill(ax, [oralabset.Startdn(a) oralabset.Stopdn(a) oralabset.Stopdn(a) oralabset.Startdn(a)], ...
                          [yl(1) yl(1) yl(2) yl(2)], oralcol, 'FaceAlpha', '0.1', 'EdgeColor', 'none');
            end
            for a = 1:size(trplset,1)
            fill(ax, [trplset.Startdn(a), trplset.Startdn(a) + 1, trplset.Startdn(a) + 1, trplset.Startdn(a)] , ...
                         [yl(1) yl(1) yl(2) yl(2)], trplcol, 'EdgeColor', 'none');
            end
            for a = 1:size(drugset,1)
                fill(ax, [drugset.Startdn(a), drugset.Startdn(a) + 1, drugset.Startdn(a) + 1, drugset.Startdn(a)], ...
                         [yl(1) yl(1) yl(2) yl(2)], drugcol, 'EdgeColor', 'none');
            end
            grid('on')
            hold off;
        end
    end
    
    if exist('f', 'var')
        savePlotInDir(f, imagefilename, subfolder);
        close(f);
    end

    toc
    
end

    
    