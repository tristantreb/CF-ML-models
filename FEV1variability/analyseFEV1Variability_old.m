% analysis of FEV1 variability with smoothing spline and movmean
%
% "Description"
% 
% Input:
% ------
% "var                  var definition (if not self-explicit)"
% "filename             what it contains (if not self-explicit)"
%
% Output:
% -------
% "variables, files, plots in the same way as input"

%% load the data

clear;
addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/smartcare/source_code;
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
study = 'BR';

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
brDrugTherapy.DateNum = addDateNum(brDrugTherapy.DrugTherapyStartDate, broffset);

% extract FEV signal

% mask for rows with FEV1 as recording type
i = ismember(brphysdata(:,5).(1), {'FEV1Recording'        });
% extract SmartCareID, DateNum, FEV
FEVdata = table2array(brphysdata(i,[1 3 8]));
clear i

% remove patient 601 and 621 (erroneous behavior) 
FEVdata = FEVdata(~ismember(FEVdata(:,1), [621 601]),:);
% patient 601 has 233 values (FEVdata(:,1) == 601) 

p_all_patients = unique(FEVdata(:,1));


%% find ideal gap

all_gaps = [];
for p_patient = p_all_patients'
    mask = FEVdata(:,1) == p_patient;
    all_gaps = [all_gaps; diff(FEVdata(mask,2))];
end
% a diff of 1 is a gap of 0 day, hence we substract 1
all_gaps = all_gaps - 1;
    
all_gaps = sort(all_gaps,'descend');
[a aa] = hist(all_gaps,100);
histogram(all_gaps(1:end-10),30)
gap_data(1,:) = 1:20; % #days gap
gap_data(2,:) = sum( all_gaps <= gap_data(1,:) );
gap_data(3,:) = gap_data(2,:)/length(all_gaps)*100;

% we allow a gap of 2 days, thus including 83% of the total measurements


%% fit the FEV curve and compute absolute values of deviation from fit

% parameters
plotfolder = '../../PlotsFEVAnalysis';
n_records_min = 50;
p_fittype=1;
p_gap = 5; 

% stable period filter parameters
p_filter = 1; % 1 for treatments and modulators, 2 for t, 3 for m only
n_prior_t = 30; % days prior to treatment start
n_post_t = 15; % days post treatment end
n_post_m = 15; % days after modulator therapy start

%% fit for all patients

i=1;
for p_smoothing = 0.1
    max_movmean_window = (p_smoothing-1)*p_gap;
%for p_smoothing = linspace(0.00001,0.001,20)
    r_all_residuals=[];
    r_all_std=[];
    r_outliers=[];
    r_patientMaxVal=[];
    r_n_used=0;
    r_entries_count=[];
    for p_patient = p_all_patients'

        % mask revealing patient data
        mask = FEVdata(:,1) == p_patient; f_datapoints = sum(mask);
        fprintf(['Patient ' num2str(p_patient) ' - ' num2str(f_datapoints) ' raw entries, ' ]);
        
        % FEVdata mask
        if p_filter ~= 0
            % mask revealing patient data during stable period
            [mask, ~] = getStableIdx(p_patient, FEVdata, ...
                ivandmeasurestable, n_prior_t, n_post_t, ...
                brDrugTherapy, n_post_m, ...
                p_filter);
            f_t_m = sum(mask);
            fprintf([num2str(f_t_m) ' stable entries, ' ]);
        end
        
        % mask revealing patient data during consecutive days (with
        % allowed gap)
        if p_fittype == 2
            [~, dates2keep] = extractConsecutiveSequences(FEVdata(mask,2)', p_gap);
            dates2remove = FEVdata(mask,2);
            dates2remove(ismember(dates2remove,dates2keep))=[];
            
            % only keep days with low gap number
            mask = not(ismember(FEVdata(:,2),dates2remove)) & mask;
        end
        fprintf([num2str(sum(mask)) ' consecutive entries (with gap).\n'])
        r_entries_count = [r_entries_count; sum(mask)];

        % check data count 
        if checkDataCount(mask, n_records_min) == 0
            r_all_std = cat(1,r_all_std, nan);
            r_patientMaxVal = [r_patientMaxVal; nan];
        else
            % count used patients
            r_n_used=r_n_used+1;

            % apply fit
            x = FEVdata(mask, 2); % date
            y = FEVdata(mask, 3); % measure
            [residuals, ~] = applyfit(x,y,p_smoothing,p_fittype);

            r_all_residuals = cat(1,r_all_residuals, residuals);
            r_all_std = cat(1,r_all_std, std(residuals));
            r_patientMaxVal = [r_patientMaxVal; max(abs(FEVdata(mask, 3))) - ...
                mean(FEVdata(mask, 3))]; % most extreme value across patients used
        end
    end

    % compute results
    a(i,1) = std(r_all_residuals);
    a(i,2) = prctile(r_all_residuals,99.5);
    a(i,3) = prctile(r_all_residuals,0.5);
    a(i,4) = prctile(r_all_residuals,97.5);
    a(i,5) = prctile(r_all_residuals,2.5);
    a(i,6) = prctile(r_all_residuals,95);
    a(i,7) = prctile(r_all_residuals,5);
    a(i,8) = sum(power(r_all_residuals,2)); % sse
    a(i,9) = p_smoothing;
    i=i+1;
end

%% segment patient types
% limit: 550 556 584
patients_perfect = p_all_patients(r_all_std < 0.063);
patients_good = p_all_patients(r_all_std >= 0.063 & r_all_std < 0.088);
patients_rough = p_all_patients(r_all_std >= 0.088 & r_all_std < 0.16);
patients_outliers = p_all_patients(r_all_std > 0.16);

%% results

figure('DefaultAxesFontSize',12,'Position', [1 1 2000 500])

subplot(3,2,2)
boxplot(r_all_std,'Orientation','horizontal');
title("Distribution of the std of residuals across patients")
xlabel('Standard deviation (L)')
grid('on')

% subplot(3,2,4)
% plotFEVModel(501, FEVdata, p_smoothing);
% grid('on')

subplot(3,2,[1,3,5])
histogram(r_all_residuals);
ylabel(['Frequency (count: ' num2str(size(r_all_residuals,1)) ...
    ' out of ' num2str(length(FEVdata(:,1))) ')'])
xlabel('Residuals (L)')
title(['Distribution of the FEV1 measures deviation from fitted curve across ' ...
    num2str(r_n_used) ' patients (out of ' ...
    num2str(length(unique(FEVdata(:,1)))) ')'], ...
    ['Std: ' num2str(std(r_all_residuals),2) ' L | ' ...
    '99% data within [' num2str(prctile(r_all_residuals,0.5),2) ', ' num2str(prctile(r_all_residuals,99.5),2) '], '...
    '95% data within [' num2str(prctile(r_all_residuals,2.5),2) ', ' num2str(prctile(r_all_residuals,97.5),2) '], '...
    '90% data within [' num2str(prctile(r_all_residuals,5),2) ', ' num2str(prctile(r_all_residuals,95),2) ']']);
grid('on')

saveas(gcf,fullfile(plotfolder,['fevModelBasedAnalysis_p' ...
    num2str(p_smoothing) '_threshold' num2str(n_records_min) ...
    '_filter' num2str(p_filter) '_fit' num2str(p_fittype) '.png']))

%% plot all curves

p_all_patients = unique(FEVdata(:,1))';
for i = 1:ceil( size(p_all_patients,2) / 10 ) % patient batch
    figure('DefaultAxesFontSize',12,'Position', [1 1 2000 1000])
    
    for j = 1:10 % 10 plots
        % get patient
        p_patient = p_all_patients(10*(i-1)+j);
        subplot(5,2,j)
        [~, days_t, days_m] = getStableIdx(p_patient, FEVdata, ...
            ivandmeasurestable, n_prior_t, n_post_t, ...
            brDrugTherapy, n_post_m, ...
            p_filter);

        plotFEVModel(p_patient, FEVdata, p_smoothing, 2, days_t, days_m); 
    end
    
    saveas(gcf,fullfile(plotfolder,['fevAllFits_p' num2str(p_smoothing) ...
        '_' num2str(p_all_patients(10*i-9)) '_fit' num2str(fittype) '.png']))
    close all
end

%% function

function plotFEVModel(patient,FEVdata,smoothing_p,threshold, days_t,days_m,fittype)
    
    p_mask = FEVdata(:,1) == patient;
    % check
    if checkDataCount(p_mask, threshold) ~= 0
        
        %fit
        x = FEVdata(p_mask, 2);
        % center data on mean
        y = FEVdata(p_mask, 3) - mean(FEVdata(p_mask, 3));
        [residuals, curve] = applyfit(x,y,smoothing_p,fittype);
        
        hold on;
        % ylim
        yl=[-0.8 0.8];
        
        % colors
        cvcol   = [0     1     0   ];
        admcol  = [0.694 0.627 0.78]; 
        ivcol   = [1     0     0   ];
        oralcol = [1     0.85  0   ];
        modcol  = [0     1     1   ];
        
        % first date
        first_date = min(x);
        
        % plot curve
        switch fittype
            case 1
                plot(curve,x,y);
            case 2
                plot(x,curve,'r',x,y,'.','MarkerEdgeColor','b');
        end
        
        xlabel('Day')
        ylabel('FEV (L)')
        ylim(yl)
        legend('off')
        title(['Patient ' num2str(patient) ' |  std: ' ...
            num2str(std(residuals),2) ' L | measures count: '...
            num2str(length(y)) ' | p: ' num2str(smoothing_p)])
        %plot(x,y,'.','b')
        
         % plot treatments influence period
         % filter days prior to start date
         days_t = days_t(days_t>first_date);
         % check nonul
         seq = extractConsecutiveSequences(days_t(days_t>first_date));
         for a = 1:size(seq,2)
            temp = seq{a};
            fill([temp(1) temp(end) temp(end) temp(1)], ...
                [yl(1) yl(1) yl(2) yl(2)], ivcol, 'FaceAlpha', '0.1', 'EdgeColor', 'none');
         end
         
         % plot modulators influence period
         seq = extractConsecutiveSequences(days_m(days_m>first_date));
         for a = 1:size(seq,2)
            temp = seq{a};
            fill([temp(1) temp(end) temp(end) temp(1)], ...
                [yl(1) yl(1) yl(2) yl(2)], modcol, 'FaceAlpha', '0.1', 'EdgeColor', 'none');
         end
        hold off;
    end
end

function is_taken = checkDataCount(data, threshold)
% boolean value to use the patient or not based on condition
    if sum(data) > threshold
        is_taken = 1;
    else
        is_taken = 0;
    end
end

function [residuals, curve] = applyfit(x, y ,p, type)
    switch type
        case 1
            [curve, goodness, output] = fit(x,y,'smoothingspline','SmoothingParam',p);
            residuals = output.residuals;
        case 3
            curve = smooth(y,p);
            residuals = y - curve;
    end
end

function out = addDateNum(date, broffset)
% create serial date referenced on broffset
    out = ceil(datenum(date) - broffset);
end