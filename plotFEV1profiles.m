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

clear;
addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/smartcare/;
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
study = 'BR';
plotfolder = getPlotFolder();

% load measures
[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[brphysdata, broffset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
% load treatments
load(fullfile(basedir, subfolder, 'BRivandmeasures_gap10.mat'));
% load CFTR modulators therapy
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brDrugTherapy');

% clean modulators tables
brDrugTherapy.DrugTherapyType = cleanDrugNamings(brDrugTherapy.DrugTherapyType);
% adds columns with serial date num
brDrugTherapy.DateNum = datenum(brDrugTherapy.DrugTherapyStartDate) - broffset;

%% extract FEV signal

% mask for rows with FEV1 as recording type
i = ismember(brphysdata(:,5).(1), {'FEV1Recording'        });
% extract SmartCareID, DateNum, FEV
FEVdata = table2array(brphysdata(i,[1 3 8]));
clear i

% parameters
n_prior_t=0;
n_post_t=0;
n_post_m=1; % min days to let a blue rectanlge appear
p_filter=1; % treatments and modulators


p_all_patients = unique(FEVdata(:,1));
patients_missing_trikafta = [835, 812, 803, 598, 525, 523, 517, 516, 501];

for patient = patients_missing_trikafta%501%p_all_patients'

    % mask revealing patient data
    mask_patient = FEVdata(:,1) == patient;

    % mask revealing patient data during stable period
    [mask, days_t, days_m] = getStableIdx(patient, FEVdata, ...
        ivandmeasurestable, n_prior_t, n_post_t, ...
        brDrugTherapy, n_post_m, ...
        p_filter);

    % check data count 
    if sum(mask_patient) >= 1 % at least one datapoint
       
        % define colors
        cvcol   = [0     1     0   ];
        admcol  = [0.694 0.627 0.78]; 
        ivcol   = [1     0     0   ];
        oralcol = [1     0.85  0   ];
        modcol  = [0     1     1   ];
        
        % plot patient
        figure('DefaultAxesFontSize',12,'Position', [1 1 2000 300])
        hold on

        % referenced on patient mean
        avg = mean(FEVdata(mask_patient, 3)); 
        ylim([-0.8 0.8]); yl=[-0.8 0.8];
        
        x = FEVdata(mask_patient, 2); % date
        y = FEVdata(mask_patient, 3); % measurements

        % plot raw measures
        plot(datetime(datestr(x+broffset)),y-avg,'.','MarkerEdgeColor','b')
        xlabel('Day')
        ylabel('FEV1 (L)')
        title(['Patient ' num2str(patient)])
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
                [yl(1) yl(1) yl(2) yl(2)], modcol, 'FaceAlpha', '1', 'EdgeColor', 'none');
         end
         
         hold off
         
         saveas(gcf,fullfile(plotfolder,['FEV1Profile_ID' num2str(patient) '.png']))
         close all
        
    end
end