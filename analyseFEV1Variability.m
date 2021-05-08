% analysis of FEV1 variability with moving average manually implemented
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
brDrugTherapy.DateNum = datenum(brDrugTherapy.DrugTherapyStartDate) - broffset;

% extract FEV signal

% mask for rows with FEV1 as recording type
i = ismember(brphysdata(:,5).(1), {'FEV1Recording'        });
% extract SmartCareID, DateNum, FEV
FEVdata = table2array(brphysdata(i,[1 3 8]));
clear i

% remove patient 601, 621 578 (erroneous behavior) 
FEVdata = FEVdata(~ismember(FEVdata(:,1), [621 601 578]),:);

p_all_patients = unique(FEVdata(:,1));

%% fit the FEV curve and compute absolute values of deviation from fit

p_all_patients = unique(FEVdata(:,1));

% parameters
plotfolder = '../../PlotsFEVAnalysis';
n_records_min = 1;
p_plot = 0;

% "moving mean"
p.window = 21; p.threshold = 7;

% stable period filter parameters
p_filter = 1; % 1 for treatments and modulators, 2 for t, 3 for m
n_prior_t = 30; % days prior to treatment start
n_post_t = 15; % days post treatment end
n_post_m = 15; % days after modulator therapy start

% %% perfect patients (from std segmentation)
perfect_patients = [501;533;538;603;639;641;603];
patients_high_dsignal_dnoise = [582, 572, 559, 558, 549, 517];
patients_high_dnoise = [802, 585, 584,523,520,503];
patients_missing_trikafta = [835, 812, 803, 598, 525, 523, 517, 516, 501];
patients_removed = [621 601 578];
% %zeros = [526,553,557,566,574,605,614,628,839,860]; % patients with 0 records after filtering

% fit for all patients

i=1;
a=[];

for w = p.window
    for t = p.threshold
        
        r_all_residuals=[];
        r_all_n_residuals=[]; 
        r_all_std=[];
        r_patientMaxVal=[];
        n_used=0;
        f_t_m = [];
        
        for patient = p_all_patients'

            % mask revealing patient data
            mask = FEVdata(:,1) == patient; mask_patient=mask ; f_datapoints = sum(mask_patient);
            fprintf('Patient %3i - %4i entries, ', patient, f_datapoints);
            
            % mask revealing patient data during stable period
            if p_filter ~= 0
                [mask, ~] = getStableIdx(patient, FEVdata, ...
                    ivandmeasurestable, n_prior_t, n_post_t, ...
                    brDrugTherapy, n_post_m, ...
                    p_filter);
                f_t_m = [f_t_m sum(mask)];
                fprintf('%4i stable entries, ',sum(mask));
            end

            % check data count 
            if checkDataCount(mask, n_records_min) == 0
                r_all_std = cat(1,r_all_std, nan);
                r_patientMaxVal = [r_patientMaxVal; nan];
                r_all_n_residuals = cat(1,r_all_n_residuals, 0);
                fprintf('   0 residuals.\n')
                
            else
                n_used=n_used+1;
                % apply fit
                x = FEVdata(mask, 2); % date
                y = FEVdata(mask, 3); % measure

                [residuals, curve] = applyMovingMean(x,y,w,t);
                fprintf('%4i residuals.\n', length(residuals(~isnan(residuals))));

                % log residuals
                r_all_residuals = cat(1,r_all_residuals, residuals);

                % log patient level data
                r_all_n_residuals = cat(1,r_all_n_residuals, sum(~isnan(residuals)));
                r_all_std = cat(1,r_all_std, std(residuals,'omitnan'));
                r_patientMaxVal = [r_patientMaxVal; max(abs(FEVdata(mask, 3))) - ...
                    mean(FEVdata(mask, 3))]; % most extreme value across patients used
                
                % plot patient
                if p_plot == 1 && sum(~isnan(residuals)) >  0 % if at least 1 residual

                    figure('DefaultAxesFontSize',12,'Position', [1 1 2000 300])
                    hold on

                    % referenced on patient mean
                    avg = mean(FEVdata(mask_patient, 3)); 
                    ylim([-0.8 0.8]);
                    % plot moving average in red
                    plot(x,curve-avg,'.','MarkerEdgeColor','r')
                    plot(x,curve-avg,'r')

                    % plot raw measures
                    plot(x,y-avg,'.','MarkerEdgeColor','b')

                    % plot data 
                    if p_filter == 1 && sum(not(mask)) > 0
                        x_na = FEVdata(mask_patient & not(mask), 2);
                        y_na = FEVdata(mask_patient & not(mask), 3);
                        plot(x_na,y_na-avg,'.','MarkerEdgeColor','g')
                    end
                    hold off

                    xlabel('Day')
                    ylabel('FEV1 (L)')
                    title(['Patient ' num2str(patient) ' |  std: ' ...
                        num2str(std(residuals,'omitnan'),3) ' L | residuals count: '...
                        num2str(sum(~isnan(residuals))) ' | window: ' num2str(w) ...
                        ' | threshold: ' num2str(t) ])
                    grid('on')
                    legend('mov avg point', 'mov avg curve','measure in stable period', 'measure in unstable period','location','best')

                    saveas(gcf,fullfile(plotfolder,['fevMovingAvg' '_patient' num2str(patient) '_w' num2str(w) '_t' num2str(t) '.png']))
                    close all
                end
            end
        end

    % compute results
    a(i,1) = 1/sqrt(length(r_all_residuals(~isnan(r_all_residuals)))); % standard error coefficient
    a(i,2) = std(r_all_residuals,'omitnan');
    a(i,3) = a(i,2)*a(i,1); % standard error of the residuals' standard deviation
    a(i,4) = prctile(r_all_residuals,99.5); % prctile treats NaNs as missing values and removes them.
    a(i,5) = prctile(r_all_residuals,0.5);
    a(i,6) = prctile(r_all_residuals,97.5);
    a(i,7) = prctile(r_all_residuals,2.5);
    a(i,8) = prctile(r_all_residuals,95);
    a(i,9) = prctile(r_all_residuals,5);
    a(i,20) = prctile(r_all_residuals,90);
    a(i,21) = prctile(r_all_residuals,10);
    a(i,10) = sum(power(r_all_residuals,2),'omitnan'); % sse
    a(i,11) = w;
    a(i,12) = t;
    a(i,13) = length(r_all_residuals(~isnan(r_all_residuals))); % n residuals computed (nan values are residuals for points that did not pass the threshold)
    a(i,14) = sum(f_t_m); %length(r_all_residuals); % sum(r_all_n_residuals); % n stable measurements = n stable residuals
    a(i,15) = length(FEVdata); % n points initially
    a(i,16) = sum(not(r_all_n_residuals==0)); % n patients with over 1 residuals, i.e. over 1 window with n_measures > threshold
    a(i,17) = n_records_min;
    a(i,18) = n_used;%length(r_all_n_residuals(~isnan(r_all_n_residuals))); % n patients over minimal data count check
    a(i,19) = length(p_all_patients); % n patients initially
    i=i+1;
    end
end

%% plot std and percntiles
x = a(:,11); %w
%x = a(:,12); %t

figure('DefaultAxesFontSize',12,'Position', [1 1 500 1000])
subplot(6,1,1)
plot(x,-a(:,13),'o--k')
title('-#residuals')
subplot(6,1,2)
plot(x,a(:,2),'o--r')
title('Std dev')
subplot(6,1,3)
plot(x(1:end-1),diff(a(:,2)),'o--r')
title('Std dev diff')
%plot(x,a(:,13),'o--r') % amount of data
subplot(6,1,4)
plot(x,a(:,4),'o--r',x,-a(:,5),'o--b')
legend('99.5','0.5','location','best')
title('Percentiles')
subplot(6,1,5)
plot(x,a(:,6),'o--r',x,-a(:,7),'o--b')
legend('97.5','2.5','location','best')
subplot(6,1,6)
plot(x,a(:,8),'o--r',x,-a(:,9),'o--b')
legend('95','5','location','best')
saveas(gcf,fullfile(plotfolder,'fevAnalysis_allpatients_w_halfw-05_filter1.png'))

%% segment patient types
% limit: 550 556 584
patients_perfect = p_all_patients(r_all_std < 0.038 & r_all_std ~= 0 & r_all_n_residuals > 100);
patients_25_50 = p_all_patients(r_all_std >= 0.045 & r_all_std < 0.0666);
patients_50_75 = p_all_patients(r_all_std >= 0.0666 & r_all_std < 0.099);
patients_outliers = p_all_patients(r_all_std > 0.18);

%% 
idx_ok = r_all_n_residuals >= 10;
kept = sum(idx_ok);
r_all_std_ok = r_all_std(idx_ok);

%% results

figure('DefaultAxesFontSize',12,'Position', [1 1 2000 500])

subplot(3,2,5)
boxplot(r_all_std,'Orientation','horizontal');
title("Boxplot of the standard deviation of patient's residuals for (21,7)", ...
    "5 out of 192 patients have 1 residual, hence a standard deviation of 0")
xlabel('Standard deviation (L)')
xticks(0:0.02:0.36)
grid('on')

% subplot(3,2,4)
% plotFEVModel(501, FEVdata, p_smoothing);
% grid('on')

subplot(3,2,2)
plot(r_all_n_residuals,'.')
title('#residuals per patient',[ num2str(sum(r_all_n_residuals==0)) ' patients with 0 residuals'])
xlabel('Patient')
ylabel('#residuals')

subplot(3,2,4)
temp = 1:31;
semilogy(temp,1./temp.*100, 'o--k')
xlabel('w')
ylabel('1/w in %')
yticks([4, 6, 8, 10, 20, 50, 100])
title("Moving mean terms' weight for window size w")
grid('on')
% put more y points

subplot(3,2,[1,3])
histogram(r_all_residuals);
ylabel(['Frequency (total: ' num2str(a(1,13)) ...
    ' from ' num2str(a(1,14)) ' filtered)'])
xlabel('Residuals (L)')
title(['Distribution of the FEV1 measures deviation from fitted curve across ' ...
    num2str(a(1,16)) ' patients (out of ' ...
    num2str(length(unique(FEVdata(:,1)))) ')'], ...
    ['Std: ' num2str(std(r_all_residuals,'omitnan'),2) ' L | ' ...
    '99% data within [' num2str(prctile(r_all_residuals,0.5),2) ', ' num2str(prctile(r_all_residuals,99.5),2) '], '...
    '95% data within [' num2str(prctile(r_all_residuals,2.5),2) ', ' num2str(prctile(r_all_residuals,97.5),2) '], '...
    '90% data within [' num2str(prctile(r_all_residuals,5),2) ', ' num2str(prctile(r_all_residuals,95),2) ']']);
grid('on')

% saveas(gcf,fullfile(plotfolder,['fevModelBasedAnalysis_movmean_w' ...
%     num2str(w) '_t' num2str(t) '_threshold' num2str(n_records_min) ...
%     '_filter' num2str(p_filter) '.png']))

%% function

function is_taken = checkDataCount(data, threshold)
% boolean value to use the patient or not based on condition
    if sum(data) >= threshold
        is_taken = 1;
    else
        is_taken = 0;
    end
end

function [residuals, curve] = applyMovingMean(x, y, w, t)
% t is inclusive
   curve = NaN(length(x),1);
   for i = 1:length(x) % to change, boundaries
       % even window: w/2 points to left, w/2-1 points to right
       % 8: [4 3]
       % odd window: w/2-0.5 points to left and right
       % 9: [4 4]
       valid_dates = x(i)-floor(w/2) : x(i)+ceil(w/2)-1;
       valid_idx = ismember(x,valid_dates);
       if sum(valid_idx) >= t
           % average values in window
           curve(i) = mean( y(valid_idx) );
       end % else value is nan by default
   end
   residuals = y - curve;
end
