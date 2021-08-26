% insights about the amount of FEV1 recordings per patient across the study

init;

% load measures
[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[brphysdata, broffset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brPatient');

%% load patient brPatient files
brPatient = [readtable('../../DataFiles/BR/PatientbrPatientFiles/PBPatientbrPatientCDF20210305.xlsx'); ...
    readtable('../../DataFiles/BR/PatientbrPatientFiles/PBPatientbrPatientPAP20210315.xlsx')];
brPatient  = brPatient(:,[1 4 6 7]);
% note: brPatient 20210315 contains 266 patients
% latest date
maxStudyDate = datenum(getLatestBreatheMeasDate,'yyyymmdd') - broffset;

%% extract FEV signal

% bluetooth recordings
% FEV1 data
% mask for rows with FEV1 as recording type
i = ismember(brphysdata(:,5).(1), {'FEV1Recording'        });
% extract SmartCareID, DateNum, FEV1
dataFEV1 = table2array(brphysdata(i,[1 3 8]));

% FEV6
i = ismember(brphysdata(:,5).(1), {'FEV6Recording'        });
dataFEV6 = table2array(brphysdata(i,[1 3 8]));

% O2_saturation
i = ismember(brphysdata(:,5).(1), {'O2SaturationRecording'        });
dataO2_saturation = table2array(brphysdata(i,[1 3 10]));

% weight
i = ismember(brphysdata(:,5).(1), {'WeightRecording'        });
dataweight = table2array(brphysdata(i,[1 3 9]));

% manual recording
% wellness
i = ismember(brphysdata(:,5).(1), {'WellnessRecording'        });
datawellness = table2array(brphysdata(i,[1 3 14]));

% cough
i = ismember(brphysdata(:,5).(1), {'CoughRecording'        });
datacough = table2array(brphysdata(i,[1 3 14]));

% temperature
i = ismember(brphysdata(:,5).(1), {'TemperatureRecording'        });
datatemperature = table2array(brphysdata(i,[1 3 12]));

% automatic recordings
% Calories
i = ismember(brphysdata(:,5).(1), {'CalorieRecording'        });
datacalories = table2array(brphysdata(i,[1 3 13]));

% resting_heart_rate
i = ismember(brphysdata(:,5).(1), {'RestingHRRecording'        });
dataresting_heart_rate = table2array(brphysdata(i,[1 3 11]));

clear i

%% specific code for fevAnalysis
p_filter=1; % activate fevAnalysis

if p_filter == 1
    % load treatments
    load(fullfile(basedir, subfolder, 'BRivandmeasures_recovery_gap10.mat'));
    % load CFTR modulators therapy
    load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'));%,'brDrugTherapy');

    % clean modulators tables
    brDrugTherapy.DrugTherapyType = cleanDrugTherapyNamings(brDrugTherapy.DrugTherapyType);
    % adds columns with serial date num
    brDrugTherapy.DateNum = datenum(brDrugTherapy.DrugTherapyStartDate) - broffset;

    n_prior_t = 30; % days prior to treatment start
    n_post_t = 15; % days post treatment end
    n_post_m = 15; % days after modulator therapy start
    
    % remove patient 601, 621 578 (erroneous behavior) 
    patients_erroneous = [221,201,178]; 
    dataFEV1 = dataFEV1(~ismember(dataFEV1(:,1), patients_erroneous),:);
end

%% build array with id, seniority and commitment

% enter the data type you want to use
type = "FEV1";
data = eval("data"+type);

patients = unique(data(:,1));
commitment = nan(length(patients),6);
mask_stable_data = zeros(length(data(:,1)),1);

for i = 1:length(patients)
    % patient id
    commitment(i,1) = patients(i);
     
    % create mask
    if p_filter == 1
        % mask revealing patient data during stable period - specific to
        % fevModel2
        [mask, days_t, days_m]  = getStableIdx(patients(i), data, ...
            ivandmeasurestable, n_prior_t, n_post_t, ...
            brDrugTherapy, n_post_m, ...
            p_filter);
        mask_stable_data = mask_stable_data | mask;
    else
        % mask revealing patient's values
        mask = data(:,1) == patients(i);
    end
    
    if sum(mask) == 0 % if no points after filter
        minDate = nan; maxDate = nan;
    else
        
        % get minDate
        % if no study start date for patient
        if isnat(brPatient.StudyDate(brPatient.ID == patients(i)))
            %  look at consent date
            if isnat(brPatient.ConsentDate(brPatient.ID == patients(i)))
                % if no consent date, take overall study start date
                minDate = 0;
            else % else take consent date
                minDate = datenum(brPatient.ConsentDate(brPatient.ID == patients(i))) - broffset;
            end
        else % else take study start date
            minDate = datenum(brPatient.StudyDate(brPatient.ID == patients(i))) - broffset;
        end
        
        % get maxDate
        if isnat(brPatient.WithdrawalDate(brPatient.ID == patients(i)))
            % if no withdrawal, take current maxDate
            maxDate = maxStudyDate;
        else % else take withdrawal date
            maxDate = datenum(brPatient.WithdrawalDate(brPatient.ID == patients(i)) - broffset);
        end
    end
    
    commitment(i,2) = minDate;
    commitment(i,3) = maxDate;
    
    % measurements' count
    commitment(i,4) = size(data(mask,:),1);
    
    % commitment: measurements' count / theoritical max in %
    if p_filter == 1 
        % remove unstable dates that are outisde from the enrollment time
        days_t = days_t(ismember(days_t, minDate:maxDate)); 
        days_m = days_m(ismember(days_m, minDate:maxDate));
        % compute number of unstable days
        n_unstable_days = length(days_t) + length(days_m);
        % divide data by theoretical max time interval
        commitment(i,6) = 100 * commitment(i,4) / (maxDate - minDate + 1 - n_unstable_days);
    else
        commitment(i,6) = 100 * commitment(i,4) / (maxDate - minDate + 1);
    end
    
    % measurement' count in % over total #measurements
    commitment(i,5) = commitment(i,4) * 100 / length(data(:,1));
end
clear i

% plot patient commitment
threshold=50;

figure('DefaultAxesFontSize',12,'Position', [1 1 500 300])
histogram(commitment(:,6),10)
hold on
histogram(commitment(commitment(:,4)>threshold,6),10,'BinLimit',[0,100])
title(sprintf('Commitment to %s recording',strrep(type,"_"," ")),  ...
    [' mean: ' num2str(mean(commitment(:,6),'omitnan'),2), '%, median: ' num2str(median(commitment(:,6),'omitnan'),2) '%, ' ...
    num2str(sum(commitment(:,4))) ' measurements, ' num2str(sum(commitment(:,4)~=0)) ' patients'])% out of ' num2str(length(unique(brphysdata.SmartCareID)))])
xlabel("Density of measurements over stable enrolment period %")
ylabel('Number of patients')
legend('all measurements', [' > ' num2str(threshold) ' measurements per patient'])
grid('on')
hold off

saveas(gcf,fullfile(plotfolder,sprintf('patientCommitment%s.png',type)))
close all
%% Create table with evolution of data count each 6 months - unused

days_span = 364/2; % 6 months window

% create 6 months date span from maxDate
day_span = minDate:days_span:maxDate; 
day_span(1)=[];
day_span(end+1)=maxDate;

t = table;

for i = 1:length(day_span)-1
    % filter data over the selected interval
    f_data = commitment(commitment(:,2) < day_span(i),:);

    % # months from study start
    t.MonthsFromStart(i) = convertCharsToStrings([num2str((i-1)*6) '-' num2str(i*6)]);
    % patient count
    t.PatientCount(i) = size(f_data,1);
    % commitment levels
    t.Commitment0to25(i) = 100*sum(0 <= f_data(:,6) & f_data(:,6) < 25)/size(f_data,1);
    t.Commitment25to50(i) = 100*sum(25 <= f_data(:,6) & f_data(:,6) < 50)/size(f_data,1);
    t.Commitment50to75(i) = 100*sum(50 <= f_data(:,6) & f_data(:,6) < 75)/size(f_data,1);
    t.Commitment75to100(i) = 100*sum(75 <= f_data(:,6) & f_data(:,6) < 100)/size(f_data,1);
end
%% find percentage of patients and measurements corresponding to a density threshold

idx = commitment(:,6)>=30; % idx patients over D% density of measurements
%idx = data(:,4)>50; % idx patients over 50 measurements count

prct_patients = 100*sum(idx)/220%length(data(:,1))*100; 

prct_measurements = 100*sum(commitment(idx,4))/sum(commitment(:,4)) % all measurements: %length(FEVdata);%

%% investigate data completeness - unused

patients_FEV1 = unique(dataFEV1(:,1));
patients_FEV6 = unique(dataFEV6(:,1));
patients_calories = unique(dataCalories(:,1));
patients_wellness = unique(dataWellness(:,1));

% patients with fev1 but no calories
patients_FEV1(not(ismember(patients_FEV1,patients_calories)))

% patients with calories but no fev1
patients_calories(not(ismember(patients_calories,patients_FEV1)))

% patients with fev1 but no wellness
patients_FEV1(not(ismember(patients_FEV1,patients_wellness)))

% patients with wellness but no calories
patients_wellness(not(ismember(patients_wellness,patients_calories)))

% patients with fev1 but no fev6
patients_FEV1(not(ismember(patients_FEV1,patients_FEV6)))