% plot FEV1 profiles with modulators and treatments
% 
% Input:
% ------
% clinical and measurements data
%
% Output:
% -------
% one plot per patient

%% load the data

init;

% load measures
[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[brphysdata, broffset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
% load treatments
load(fullfile(basedir, subfolder, 'BRivandmeasures_gap10.mat'));
% load CFTR modulators therapy
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brDrugTherapy');

% clean modulators tables
brDrugTherapy.DrugTherapyType = cleanDrugTherapyNamings(brDrugTherapy.DrugTherapyType);
% adds columns with serial date num
brDrugTherapy.DateNum = datenum(brDrugTherapy.DrugTherapyStartDate) - broffset;

%% extract signals

% mask for rows with FEV1 as recording type
i = ismember(brphysdata(:,5).(1), {'FEV1Recording'        });
% extract SmartCareID, DateNum, FEV1
dataFEV1 = table2array(brphysdata(i,[1 3 8]));

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

% cold or flue
i = ismember(brphysdata(:,5).(1), {'HasColdOrFluRecording'});
dataresting_cold_flue = table2array(brphysdata(i,[1 3 16]));
clear i

%% 
% parameters
n_prior_t=0;
n_post_t=0;
n_post_m=1; % take one day for the modulators
p_filter=1; % treatments and modulators


p_all_patients = 501%unique(FEVdata(:,1));
patients_missing_trikafta = [835, 812, 803, 598, 525, 523, 517, 516, 501];

for patient = p_all_patients'

    % mask revealing patient data
    mask_patient = dataFEV1(:,1) == patient;

    % mask revealing patient data during stable period
    [mask, days_t, days_m] = getStableIdx(patient, dataFEV1, ...
        ivandmeasurestable, n_prior_t, n_post_t, ...
        brDrugTherapy, n_post_m, ...
        p_filter);

    % check data count 
    if sum(mask_patient) >= 1 % at least one datapoint
       
        % define colors
        cvcol   = [0.94  0.52  0.15];
        admcol  = [0.694 0.627 0.78]; 
        ivcol   = [1     0     0   ];
        oralcol = [1     0.85  0   ];
        trplcol = [0     1     0   ];
        drugcol = [0     0.8   0.6 ];
        
        % plot patient
        figure('DefaultAxesFontSize',12,'Position', [1 1 2000 300])
        
        subplot(8,1,1)
        
        % referenced on patient mean
        avg = mean(dataFEV1(mask_patient, 3)); 
        ylim([-0.8 0.8]); yl=[-0.8 0.8];
        
        x = dataFEV1(mask_patient, 2); % date
        y = dataFEV1(mask_patient, 3); % measurements

        % plot raw measures
        plot(datetime(datestr(x+broffset)),y-avg,'.','MarkerEdgeColor','b')
        xlabel('Day')
        ylabel('FEV1 (L)')
        title(sprintf('Patient %i, mean FEV1: %1.1f L', patient, avg))
        grid('on')

        % first patient's measurement date
        first_date = min(x);

         % plot IV & oral influence period
         seq = extractConsecutiveSequences(days_t(days_t>first_date),0);
         for a = 1:size(seq,2)
            temp = datetime(datestr(seq{a}+broffset));
            fill([temp(1) temp(end) temp(end) temp(1)], ...
                [yl(1) yl(1) yl(2) yl(2)], ivcol, 'FaceAlpha', '0.1', 'EdgeColor', 'none');
         end

         % plot modulators influence period
         seq = extractConsecutiveSequences(days_m(days_m>first_date),0);
         for a = 1:size(seq,2)
            temp = datetime(datestr(seq{a}+broffset));
            fill([temp(1) temp(end) temp(end) temp(1)], ...
                [yl(1) yl(1) yl(2) yl(2)], trplcol, 'FaceAlpha', '1', 'EdgeColor', 'none');
         end
         
         %saveas(gcf,fullfile(plotfolder,['FEV1Profile_ID' num2str(patient) '.png']))
         %close all
        
    end
end

%function plot