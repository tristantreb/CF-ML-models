% analysis of FEV1 variability with moving average manually implemented
%
% Model based statistical estimation of the variability in FEV1 measurements.

%% load the data

init;
plotfolder = '../../PlotsFEVAnalysis';

% load measures
[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[brphysdata, broffset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
% load treatments
load(fullfile(basedir, subfolder, 'BRivandmeasures_recovery_gap10.mat'));
% load CFTR modulators therapy
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brDrugTherapy','brPatient');

%% clean modulators tables
%brDrugTherapy.DrugTherapyType = cleanDrugTherapyNamings(brDrugTherapy.DrugTherapyType);
% adds columns with serial date num
brDrugTherapy.DateNum = datenum(brDrugTherapy.DrugTherapyStartDate) - broffset;
brDrugTherapy.StopDateNum = datenum(brDrugTherapy.DrugTherapyStopDate) - broffset;
% remove patients with negative dates
neg_dates = size(brDrugTherapy,1);
brDrugTherapy(brDrugTherapy.DateNum < 0,:)=[];
neg_dates = neg_dates - size(brDrugTherapy,1);

%% extract FEV signal

% mask for rows with FEV1 as recording type
FEVdata = getMeasureTable(brphysdata,'FEV1Recording' ,'FEV');
FEVdata = removevars(FEVdata, {'UserName', 'CaptureType', 'Date_TimeRecorded', 'ScaledDateNum'});

% remove patient 601, 621 578 (erroneous behavior) 
%patients_erroneous = [621,601,578]; old ID
patients_erroneous = [221,201,178]; 
patients_outliers = [115;141;178;196;201;235;275;279]; %6 more than for erroneous patients
FEVdata = FEVdata(~ismember(FEVdata.ID, patients_erroneous),:);
FEVdata = FEVdata(~ismember(FEVdata.ID, patients_outliers),:);

% this code was done with FEVdata as an array, transforming
FEVdata = table2array(FEVdata);

%% fit the FEV curve and compute absolute values of deviation from fit

% parameters
p_plot = 0;
n_records_min = 1;
p.window = 21; p.threshold = 7;
% filter stable period
p_filter = 1; % 1 to remove days during stable period, 2/3 to base it on treatments/modulators only

% stable period filter parameters
n_prior_t = 30; % days prior to treatment start
n_post_t = 15; % days post treatment end
n_post_m = 15; % days after modulator therapy start

% %% perfect patients (from std segmentation)
% perfect_patients = [501;533;538;603;639;641;603];
% patients_high_dsignal_dnoise = [582, 572, 559, 558, 549, 517];
% patients_high_dnoise = [802, 585, 584,523,520,503];
% patients_missing_trikafta = [835, 812, 803, 598, 525, 523, 517, 516, 501];
% patients_removed = [621 601 578];
% %zeros = [526,553,557,566,574,605,614,628,839,860]; % patients with 0 records after filtering

p_processed_patients = unique(FEVdata(:,1));

rPatientStats=initPatientStatTable(length(p_processed_patients));
rPatientStats.p_processed_patients = p_processed_patients;

% fit for all patients
i=1;
r_statistics=[];

for w = p.window
    for t = p.threshold
        
        n_used=0;
        f_t_m = [];
        mask_stable_all = zeros(length(FEVdata(:,1)),1); % mask with all stable measures
        r_all_residuals = [];
        
        % patients that had triple therapy
        n_1_trpl = 0;
        r_all_1_prior_trpl = []; 
        r_all_1_post_trpl = [];
        
        % patients that had symkevi and then triple therapy
        % symkevi
        n_2_smkv_and_trpl = 0;
        r_all_2_none = []; % no symkevi, no triple therapy
        r_all_2_smkv = []; % symkevi, no triple therapy
        r_all_2_trpl = []; % triple therapy
        
        % patients that had symkevi and no triple therapy
        n_3_smkv = 0;
        r_all_3_prior_smkv = [];
        r_all_3_post_smkv = [];
        
        n_residual_no_trpl = 0;
        
        for patient = rPatientStats.p_processed_patients'

            patient_idx = p_processed_patients == patient;
            % mask revealing patient data
            mask = FEVdata(:,1) == patient; mask_patient=mask;
            fprintf('Patient %3i - %4i entries', patient, sum(mask_patient));
            
            % mask revealing patient data during stable period
            if p_filter ~= 0
                [mask, ~] = getStableIdx(patient, FEVdata, ...
                    ivandmeasurestable, n_prior_t, n_post_t, ...
                    brDrugTherapy, n_post_m, ...
                    p_filter);
                f_t_m = [f_t_m sum(mask)];
                fprintf(', %4i stable entries',sum(mask));
                
                % create mask with all stable data over all patients
                mask_stable_all = mask_stable_all | mask;
            end

            % check to use the patient or not
            if sum(mask) == 0 % check non null post filtering data count

                rPatientStats.all_std(patient_idx) = nan;
                rPatientStats.patientMaxVal(patient_idx) = nan;
                rPatientStats.all_n_residuals(patient_idx) = 0;
                % used for the patient level analysis of the effect of
                % triple therapy
                rPatientStats.n_residuals_prior_tripleT(patient_idx) = nan; 
                rPatientStats.n_residuals_post_tripleT(patient_idx) = nan;
                rPatientStats.std_residuals_prior_tripleT(patient_idx) = nan; 
                rPatientStats.std_residuals_post_tripleT(patient_idx) = nan;
                fprintf(',    0 residuals.\n')
                
            else % patient is used
                n_used=n_used+1;
                % apply fit
                x = FEVdata(mask, 2); % date
                y = FEVdata(mask, 3); % measure

                [residuals, curve] = applyMovingMean(x,y,w,t);
                fprintf(', %4i residuals', length(residuals(~isnan(residuals))));

                % log residuals
                r_all_residuals = cat(1,r_all_residuals, residuals);

                % log patient level data
                rPatientStats.all_n_residuals(patient_idx) = sum(~isnan(residuals));
                rPatientStats.all_std(patient_idx) = std(residuals,'omitnan');
                rPatientStats.patientMaxVal(patient_idx) = max(abs(FEVdata(mask, 3))) - ...
                    mean(FEVdata(mask, 3)); % most extreme value across patients used
                
                
                if isempty(rPatientStats.all_n_residuals(rPatientStats.p_processed_patients == patient))% if no residuals
                    rPatientStats.n_residuals_prior_tripleT(patient_idx) = nan; 
                    rPatientStats.n_residuals_post_tripleT(patient_idx) = nan; 
                    rPatientStats.std_residuals_prior_tripleT(patient_idx) = nan;
                    rPatientStats.std_residuals_post_tripleT(patient_idx) = nan;
                else
                    
                    %%% 1 %%% split residuals prior/post triple therapy start
                    % avoid residuals after triple therapy stop if
                    % applicable

                    % get triple therapy start and stop dates where applicable
                    triple_therapy_start = brDrugTherapy.DateNum(brDrugTherapy.ID == patient & ismember(brDrugTherapy.DrugTherapyType, 'Kaftrio/Trikafta/TripleTherapy'));
                    triple_therapy_stop = brDrugTherapy.StopDateNum(brDrugTherapy.ID == patient & ismember(brDrugTherapy.DrugTherapyType, 'Kaftrio/Trikafta/TripleTherapy'));
                    if ~isempty(triple_therapy_start)
                        fprintf(', with triple therapy');
                        n_1_trpl = n_1_trpl+1;
                        if length(triple_therapy_start) > 1
                            triple_therapy_start = triple_therapy_start(1);
                        end
                        
                        % segment residuals prior/post triple therapy start
                        % prior
                        temp = residuals(x < triple_therapy_start);
                        temp = temp(~isnan(temp));
                        rPatientStats.n_residuals_prior_tripleT(patient_idx) = length(temp);
                        rPatientStats.std_residuals_prior_tripleT(patient_idx) = std(temp);
                        r_all_1_prior_trpl = cat(1,r_all_1_prior_trpl, temp);
                        
                        % post
                        if ~isempty(triple_therapy_stop)
                            temp = residuals(x >= triple_therapy_start);
                        else % avoiding measurementss after triple therapy stop
                            temp = residuals(x >= triple_therapy_start & x < triple_therapy_stop);
                        end
                        temp = temp(~isnan(temp));
                        rPatientStats.n_residuals_post_tripleT(patient_idx) = length(temp);
                        rPatientStats.std_residuals_post_tripleT(patient_idx) = std(temp); 
                        r_all_1_post_trpl = cat(1,r_all_1_post_trpl, temp); 
                        clear temp;
                        
%                         % hypothesis test that the two std are the same
%                         [h_rejected, p_value] = vartest2(r_all_prior_trpl(~isnan(r_all_prior_trpl)),r_all_post_trpl(~isnan(r_all_post_trpl)),'Tail','right');
%                         
%                         rPatientStats.h_rejected(patient_idx) = h_rejected;
%                         rPatientStats.h_p(patient_idx) = p_value; 
                        

                        %%% 2 %%% segment residuals 1) before symkevi 2)
                        %%% during symkevi, 3) during triple therapy
                            % avoiding residuals for patients that came back on
                            % symkevi - only 3 patients at the moment
                            % (18.08.2021), will not impact results
                            
                        symkevi_start = brDrugTherapy.DateNum(brDrugTherapy.ID == patient & ismember(brDrugTherapy.DrugTherapyType, 'Symkevi'));
                        if ~isempty(symkevi_start)
                            if symkevi_start(1) < triple_therapy_start
                                fprintf(' and symkevi');
                                
                                symkevi_start = symkevi_start(1);
                                n_2_smkv_and_trpl = n_2_smkv_and_trpl+1;

                                % prior symkevi
                                temp = residuals(x < triple_therapy_start & x < symkevi_start);
                                r_all_2_none = cat(1,r_all_2_none, temp(~isnan(temp)));

                                % during symkevi
                                % note: if patients started triple therapy before
                                % symkevi, we consider the symkevi data as
                                % triple therapy
                                temp = residuals(x < triple_therapy_start & x >= symkevi_start);
                                r_all_2_smkv = cat(1,r_all_2_smkv, temp(~isnan(temp)));

                                % during triple therapy
                                if ~isempty(triple_therapy_stop)
                                    temp = residuals(x >= triple_therapy_start);
                                else % avoiding measurementss after triple therapy stop
                                    temp = residuals(x >= triple_therapy_start & x < triple_therapy_stop);
                                end
                                r_all_2_trpl = cat(1,r_all_2_trpl, temp(~isnan(temp)));
                                clear temp;
                            end
                        end
                    else 
                        n_residual_no_trpl = n_residual_no_trpl + 1;
                        % update values only prior to trikafta
%                         rPatientStats.n_residuals_prior_tripleT(patient_idx) = length(~isnan(residuals)); 
%                         rPatientStats.std_residuals_prior_tripleT(patient_idx) = std(residuals,'omitnan');
%                         
%                         rPatientStats.n_residuals_post_tripleT(patient_idx) = 0; 
%                         rPatientStats.std_residuals_post_tripleT(patient_idx) = nan;
%                     
                        %%% 3 %%% semgent residuals prior/post symkevi start
                        % get symkevi start date where applicable
                        % those patients did not have triple therapy
                            % avoiding residuals for patients that came back on
                            % symkevi - only 3 patients at the moment
                            % (18.08.2021), will not impact results
                        
                        symkevi_start = brDrugTherapy.DateNum(brDrugTherapy.ID == patient & ismember(brDrugTherapy.DrugTherapyType, 'Symkevi'));
                        if ~isempty(symkevi_start)
                            fprintf(', with symkevi');
                            n_3_smkv = n_3_smkv+1;

                            % find residuals prior/post triple therapy start
                            % prior
                            temp = residuals(x < symkevi_start);
                            r_all_3_prior_smkv = cat(1,r_all_3_prior_smkv, temp(~isnan(temp)));

                            % post
                            temp = residuals(x >= symkevi_start);
                            r_all_3_post_smkv = cat(1,r_all_3_post_smkv, temp(~isnan(temp))); 
                            clear temp;
                        else
%                             fprintf(', no drug therapy');
%                             temp = residuals;
%                             r_all_3_prior_smkv = cat(1,r_all_3_prior_smkv, temp(~isnan(temp)));
                        end
                    end
                end
                fprintf('.\n');
                
                % plot patient
                if p_plot == 1 && sum(~isnan(residuals)) >  0 % if at least 1 residual

                    figure('DefaultAxesFontSize',12,'Position', [1 1 2000 300])
                    hold on

                    % referenced on patient mean
                    avg = mean(FEVdata(mask_patient, 3)); 
                    ylim([-0.8 0.8]);
                    % plot moving average in red
                    plot(x,curve-avg,'.','MarkerEdgeColor','r','MarkerSize',20)
                    %plot(x,curve-avg,'r')

                    % plot raw measures
                    plot(x,y-avg,'.','MarkerEdgeColor','b','MarkerSize',10)

                    % plot data 
                    if p_filter == 1 && sum(not(mask)) > 0
                        x_na = FEVdata(mask_patient & not(mask), 2);
                        y_na = FEVdata(mask_patient & not(mask), 3);
                        plot(x_na,y_na-avg,'.','MarkerEdgeColor',[0.8 0.8 0.9],'MarkerSize',10 )
                    end
                    hold off

                    xlabel('Day')
                    ylabel('Normalised FEV1')
                    title(['Patient ' num2str(patient) ' |  std: ' ...
                        num2str(std(residuals,'omitnan'),3) ' L | residuals count: '...
                        num2str(sum(~isnan(residuals))) ' | window: ' num2str(w) ...
                        ' | threshold: ' num2str(t) ])
                    grid('on')
                    legend('residuals','data in stable period', 'data in unstable period','location','best')

                    saveas(gcf,fullfile(plotfolder,['fevMovingAvg' '_patient' num2str(patient) '_w' num2str(w) '_t' num2str(t) '.png']))
                    close all
                end
            end
        end

    % compute results
    r_statistics(i,2) = std(r_all_residuals,'omitnan'); n = sum(rPatientStats.all_n_residuals); alpha = 0.05;
    r_statistics(i,3) = sqrt((n-1)./chi2inv(1-alpha/2,n-1))*r_statistics(i,2); % lower bound of the CI for the std
    r_statistics(i,4) = sqrt((n-1)./chi2inv(alpha/2,n-1))*r_statistics(i,2); % upper bound of the CI
    r_statistics(i,5) = prctile(r_all_residuals,99.5); % prctile treats NaNs as missing values and removes them.
    r_statistics(i,6) = prctile(r_all_residuals,0.5);
    r_statistics(i,7) = prctile(r_all_residuals,97.5);
    r_statistics(i,8) = prctile(r_all_residuals,2.5);
    r_statistics(i,9) = prctile(r_all_residuals,95);
    r_statistics(i,10) = prctile(r_all_residuals,5);
    r_statistics(i,11) = prctile(r_all_residuals,90);
    r_statistics(i,12) = prctile(r_all_residuals,10);
    r_statistics(i,13) = sum(power(r_all_residuals,2),'omitnan'); % sse
    r_statistics(i,14) = w;
    r_statistics(i,15) = t;
    r_statistics(i,16) = sum(rPatientStats.all_n_residuals); % alsolength(r_all_residuals(~isnan(r_all_residuals))); % n residuals computed (nan values are residuals for points that did not pass the threshold)
    r_statistics(i,17) = sum(f_t_m);  % n residuals during stable period
    r_statistics(i,18) = length(FEVdata); % n points initially
    r_statistics(i,19) = sum(rPatientStats.all_n_residuals > 0); % n patients with over 1 residuals, i.e. over 1 window with n_measures > threshold
    r_statistics(i,20) = n_records_min;
    r_statistics(i,21) = n_used;%length(r_all_n_residuals(~isnan(r_all_n_residuals))); % n patients over minimal data count check (stable)
    r_statistics(i,22) = length(rPatientStats.p_processed_patients); % n patients initially
    i=i+1;
    end
end

rPatientStats.s1s2 = rPatientStats.std_residuals_prior_tripleT.^2 ./ rPatientStats.std_residuals_post_tripleT.^2;

%% plot std and percntiles
x = r_statistics(:,11); %w
%x = a(:,12); %t

figure('DefaultAxesFontSize',12,'Position', [1 1 500 1000])
subplot(6,1,1)
plot(x,-r_statistics(:,13),'o--k')
title('-#residuals')
subplot(6,1,2)
plot(x,r_statistics(:,2),'o--r')
title('Std dev')
subplot(6,1,3)
plot(x(1:end-1),change(r_statistics(:,2)),'o--r')
title('Std dev diff')
%plot(x,a(:,13),'o--r') % amount of data
subplot(6,1,4)
plot(x,r_statistics(:,4),'o--r',x,-r_statistics(:,5),'o--b')
legend('99.5','0.5','location','best')
title('Percentiles')
subplot(6,1,5)
plot(x,r_statistics(:,6),'o--r',x,-r_statistics(:,7),'o--b')
legend('97.5','2.5','location','best')
subplot(6,1,6)
plot(x,r_statistics(:,8),'o--r',x,-r_statistics(:,9),'o--b')
legend('95','5','location','best')
saveas(gcf,fullfile(plotfolder,'fevAnalysis_allpatients_w_halfw-05_filter1.png'))

%% segment patient types depending on patient's variability
% limit: 550 556 584
patients_perfect = p_processed_patients(rPatientStats.all_std < 0.038 & rPatientStats.all_std ~= 0 & rPatientStats.all_n_residuals > 100);
patients_25_50 = p_processed_patients(rPatientStats.all_std >= 0.045 & rPatientStats.all_std < 0.0666);
patients_50_75 = p_processed_patients(rPatientStats.all_std >= 0.0666 & rPatientStats.all_std < 0.099);
patients_outliers = p_processed_patients(rPatientStats.all_std > 0.19);

%% results

figure('DefaultAxesFontSize',12,'Position', [1 1 2000 500])

subplot(3,2,5)
boxplot(rPatientStats.all_std,'Orientation','horizontal');
title("Boxplot of the standard deviation of patient's residuals for (21,7)", ...
    sprintf("%i out of %i patients have 1 residual, hence a standard deviation of 0", sum(rPatientStats.all_std==0), size(rPatientStats,1) ));
xlabel('Standard deviation (L)')
xticks(0:0.02:0.36)
grid('on')

% subplot(3,2,4)
% plotFEVModel(501, FEVdata, p_smoothing);
% grid('on')

subplot(3,2,2)
plot(rPatientStats.all_n_residuals,'.')
title('#residuals per patient',[ num2str(sum(rPatientStats.all_n_residuals==0)) ' patients with 0 residuals'])
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
ylabel(['Frequency (total: ' num2str(r_statistics(1,16)) ...
    ' out of ' num2str(r_statistics(1,17)) ' stable)'])
xlabel('Residuals (L)')
title(['Distribution of the FEV1 measures deviation from fitted curve across ' ...
    num2str(r_statistics(1,19)) ' patients (out of ' ...
    num2str(r_statistics(1,22)) ')'], ...
    ['Std: ' num2str(std(r_all_residuals,'omitnan'),2) ' L | ' ...
    '99% data within [' num2str(prctile(r_all_residuals,0.5),2) ', ' num2str(prctile(r_all_residuals,99.5),2) '], '...
    '95% data within [' num2str(prctile(r_all_residuals,2.5),2) ', ' num2str(prctile(r_all_residuals,97.5),2) '], '...
    '90% data within [' num2str(prctile(r_all_residuals,5),2) ', ' num2str(prctile(r_all_residuals,95),2) ']']);
grid('on')

saveas(gcf,fullfile(plotfolder,sprintf('fevModelBasedAnalysis_movmean_w%i_t%i_filter_%i.png', ...
    w, t, p_filter)))

%% FURTHER ANALYSES

%% Distribution of the residuals
% check if the assumption of normal distribution is acceptable

figure
probplot('normal',r_all_residuals)
grid('on')
figure
probplot('extreme value',r_all_residuals)
grid('on')

% better fit for a normal distribution than an extreme value distribution

% the distribution has heavy tails on both sides: probably close to a
% Cauchy distribution

% quite close to the normal distribution in the boundaries [0.05,0.95]
% which is sufficient to assume normal distribution for hypothesis testing 
% with a significance threshold of 0.05

%% plot correlation between lung health and variability

% adds a columm to FEVdata_table with FEV1 in % of predicted value
    
% Input:
% ------
% FEVdata_table        contains ID, DateNum, FEV1
% brPatient            contains ID, CalcPredictedFEV1 (this is the estimated
%   FEV1 using a function from Andres - we don't trust the PredictedFEV1  
%   that was manually entered, with errors, by the hospitals nurses)
        
% take all values again
i = ismember(brphysdata(:,5).(1), {'FEV1Recording'        });
FEV1info = brphysdata(i,[1 8]); FEV1info.Properties.VariableNames = cell({'ID','FEV1'}); clear i;

% take stable values only
% FEV1info = table(FEVdata(mask_stable_data,1),FEVdata(mask_stable_data,3)); FEV1info.Properties.VariableNames={'ID','FEV1'};

% calculate true value based on mean
func = @(x) mean(x);
FEV1info = varfun(func,FEV1info,'GroupingVariables','ID');

% adds predicted FEV1
FEV1info = outerjoin(FEV1info, brPatient, 'Type', 'Left', 'Keys', 'ID', 'RightVariables', {'CalcPredictedFEV1','Age','Height'});

% computes % predicted
FEV1info.PercentagePredicted = FEV1info.Fun_FEV1 ./ FEV1info.CalcPredictedFEV1 * 100;

% outerjoin to add patients' std and #residuals
FEV1info = outerjoin(rPatientStats, FEV1info, 'Type', 'Left', 'LeftKeys','p_processed_patients','RightKeys', 'ID', 'LeftVariables', {'all_std','all_n_residuals'});

% filter 
FEV1info = FEV1info( FEV1info.all_std <= 0.3 & FEV1info.all_n_residuals >= 30 ,:);

% variability
scatterhist(FEV1info.all_std,FEV1info.PercentagePredicted,10)
xlabel('\sigma_{residuals} (L)')
ylabel('FEV1 %')
[a,b] = corr(FEV1info.all_std,FEV1info.PercentagePredicted,'Rows','complete')
title(sprintf('Predicted FEV1 %% against patient-specific variability (r = %2.3f, p-value = %2.3f)',a,b),...
    sprintf('%i patients, %i outliers removed (\\sigma_{residuals} > 0.3 or #residuals < 30)', size(FEV1info,1), sum(not(rPatientStats.all_n_residuals==0))-size(FEV1info,1)))
clear a; clear b;

saveas(gcf,fullfile(plotfolder,'fevModelBasedAnalysis_variabilityvsFEV1.png'))
%close all

%% Effect of triple therapy at a patient level

% There is 95% probability that the CI contains the true standard deviation, 
% and if it does contain it, with 200 points, the true std value lies
% within 90-110% of the estimated std

% Hence, true variability will lie within 0.9^4-1.1^4 = 65.6-1.46% of the estimated std
% We define the variability as 4 sigma deviations, thus keeping 93% percent of the datapoints,
% according to the Chebyshev's inequality. Noise over 4 sigma is considered as outliers.

% Other option: test the difference between sigma_prior and sigma_post.
% However, those sigmas are estimators. Hence, I would have to take the
% extreme case: the boundary of the CI to be maximise the significance

% Other option: add patient-specific samples together.
% compare residuals before and after. I should have enough samples to have
% narrow CI.

varChange = rPatientStats(rPatientStats.n_residuals_prior_tripleT > 200 & rPatientStats.n_residuals_post_tripleT > 200,:);
a=varChange.std_residuals_prior_tripleT;
b=varChange.std_residuals_post_tripleT;
change=varChange.std_residuals_prior_tripleT - varChange.std_residuals_post_tripleT;
varChange.relative_change_percentage = 100*(varChange.std_residuals_prior_tripleT - varChange.std_residuals_post_tripleT) ./ varChange.std_residuals_prior_tripleT;

% histogram(varChange.relative_change_percentage,30);
% % conditions: should have start triple therapy
% % should have residuals
% % should have residuals before and after
% title(sprintf('Variability change after triple therapy start'),...
%     sprintf('%i patients have at least 2 residuals before/after (min. for a nonzero standard deviation)', sum(idx_diff)));
% xlabel('(\sigma_{residuals}^{before} - \sigma_{residuals}^{after}) / \sigma_{residuals}^{before} (L)');
% %xticks(-0.16:0.02:0.14)
% ylabel('#patients');
% grid('on')

temp_conditions = outerjoin(rPatientStats,brDrugTherapy,'Type','Left','LeftKeys','p_processed_patients','RightKeys','ID','LeftVariables','all_n_residuals','RightVariables','DrugTherapyType');
temp_conditions = temp_conditions(temp_conditions.DrugTherapyType == "Triple Therapy" & temp_conditions.all_n_residuals > 0,:);
fprintf('Variability change after triple therapy start - %i patients started triple therapy & have at least one residual, only %i patients have more than 200 residuals.\n',...
    size(temp_conditions,1),size(varChange,1)); clear temp_conditions;

saveas(gcf,fullfile(plotfolder,'fevModelBasedAnalysis_variabilityChangeAfterTripleTherapyStart.png'))
%close all

%% 1) Effect of triple therapy at a study level
% remove patients with standard deviation outliers 

distname = 'tLocationScale';

% most conservative scenario
% 1. equilibrium of data prior/after (residuals of patient that did not get
% triple therapy is not taken into account, otherwise balance largely in
% favor of prior triple therapy)
% 2. some patients got triple therapy but it was not referenced in our data,
% we simply don't use those patient by strictly looking at the ones that
% got triple therapy
% 3. symkevi data is included into the "prior" data. Since symkevi has
% potentially also had an effect, 

% 1. null hypothesis: the two std dev are the same. std_prior^2 - std_post^2 = 0
s1_prior_trpl = std(r_all_1_prior_trpl);
s1_post_trpl = std(r_all_1_post_trpl);
% relative change in s
% s1s2 = (s1-s2)/s1;

% 2. alternative hypothesis std_prior^2 - std_post^2 > 0
% assumptions: 
% . sample comes from a normally distributed population
% . measurements are independently and identically distributed, hence so
% are residuals

% test chi square distribution
% chirnd = chi2rnd(4, size(r_all_1_prior_trpl));

% data normality
figure('DefaultAxesFontSize',12,'Position', [1 1 800 300])
subplot(1,2,1)
% probplot('normal',r_all_1_prior_trpl)
pd = fitdist(r_all_1_prior_trpl,distname)
qqplot(r_all_1_prior_trpl, pd)
ylabel(sprintf('Quantiles of residuals prior Trikafta (%d)',length(r_all_1_prior_trpl) ))
title('')
grid('on')
subplot(1,2,2)
% probplot('normal',r_all_1_post_trpl)
pd = fitdist(r_all_1_post_trpl,distname)
qqplot(r_all_1_post_trpl, pd)
title('')
ylabel(sprintf('Quantiles of residuals during Trikafta (%d)',length(r_all_1_post_trpl) ))
grid('on')
saveas(gcf,fullfile(plotfolder,'fevModelBasedAnalysis_qqplot_prior_during_tripletherapy.png'))

% symmetric, heavy-tailed distribution -> Student, Cauchy

% 3. test statistic
% F-test: extremely sensitive to departures from normality for small alpha
% levels (lower than 0.05)
% Bartlett: Chi-Square, sensitive to departues from normality.
% Levene (mean) - Brown-Forsythe (median): mean is best when underlying 
% data follow Cauchy distribution, median is best for Chi-Square distribution

alpha = 0.05;

% F-Test
% F = s1^2/s2^2;
[F_h_rejected_1, F_p_1, F_ci_1, F_stats_1] = vartest2(r_all_1_prior_trpl,r_all_1_post_trpl,'Tail','right','Alpha',alpha);

% Levene test
[L_h_rejected_1, L_p_1, L_stat_1, L_critical_1] = levenetest(r_all_1_prior_trpl,r_all_1_post_trpl,alpha);

%% 2) Effect of symkevi and triple therapy at a study level

s2_smkv = std(r_all_2_smkv);
s2_trpl = std(r_all_2_trpl);
s2_none = std(r_all_2_none);

% data normality
figure('DefaultAxesFontSize',12,'Position', [1 1 1000 300])
subplot(1,3,1)
pd = fitdist(r_all_2_none,distname)
qqplot(r_all_2_none, pd)
title('')
ylabel(sprintf('Quantiles of residuals prior Symkevi (%d)',length(r_all_2_none) ))
xlabel('Quantiles of tlocationscale Distribution (\nu = 3)')
grid('on')
subplot(1,3,2)
pd = fitdist(r_all_2_smkv,distname)
qqplot(r_all_2_smkv, pd)
title('')
ylabel(sprintf('Quantiles of residuals during Symkevi (%d)',length(r_all_2_smkv) ))
xlabel('Quantiles of tlocationscale Distribution (\nu = 4)')
grid('on')
subplot(1,3,3)
pd = fitdist(r_all_2_trpl,distname)
qqplot(r_all_2_trpl, pd)
title('')
ylabel(sprintf('Quantiles of residuals during Trikafta (%d)',length(r_all_2_trpl) ))
xlabel('Quantiles of tlocationscale Distribution (\nu = 4)')
grid('on')
saveas(gcf,fullfile(plotfolder,'fevModelBasedAnalysis_qqplot_none_smkv_trpl.png'))

% symmetric, heavy-tailed distribution -> Student, Cauchy
% note: prior to Symkevi and post Symkevi (during Triple Therapy), higher
% deviation from normality around the 5th and 95th percentiles

% F-tests
[F_h_rejected_smkv_none, F_p_smkv_none, F_ci_smkv_none, F_stats_smkv_none] = vartest2(r_all_2_none,r_all_2_smkv,'Tail','right','Alpha',alpha);
[F_h_rejected_trpl_smkv, F_p_trpl_smkv, F_ci_trpl_smkv, F_stats_trpl_smkv] = vartest2(r_all_2_smkv,r_all_2_trpl,'Tail','right','Alpha',alpha);
[F_h_rejected_trpl_none, F_p_trpl_none, F_ci_trpl_none, F_stats_trpl_none] = vartest2(r_all_2_none,r_all_2_trpl,'Tail','right','Alpha',alpha,'Alpha',alpha);

% Levene's tests
[L_h_rejected_smkv_none, L_p_smkv_none, L_stat_smkv_none, L_critical_smkv_none] = levenetest(r_all_2_none,r_all_2_smkv,alpha);
[L_h_rejected_trpl_smkv, L_p_trpl_smkv, L_stat_trpl_smkv, L_critical_trpl_smkv] = levenetest(r_all_2_smkv,r_all_2_trpl,alpha);
[L_h_rejected_trpl_none, L_p_trpl_none, L_stat_trpl_none, L_critical_trpl_none] = levenetest(r_all_2_none,r_all_2_trpl,alpha);

%% 3) Effect of symkevi at a study level

% as we cannot reject the homoscedasticity hypotheses on the effects of the
% 1. transition from nothing (or Ivacaftor, or Okrambi) to Symkevi
% 2. transition from Symkevi to Triple Therapy 
% it may be because of a lack of data (only 61 patients got Symkevi and
% Triple Therapy)
% hence we have a look at Symkevi only and Symkevi & Triple Therapy
% patients

% merge patient that had symkevi and not trikafta with patient that had
% symkevi and trikafta - smae for no therapy
r_all_3_prior_smkv_full = [r_all_3_prior_smkv;r_all_2_none];
r_all_3_post_smkv_full = [r_all_3_post_smkv;r_all_2_smkv];

s3_prior_smkv = std(r_all_3_prior_smkv_full);
s3_prior_smkv2 = std(r_all_3_prior_smkv);
s3_post_smkv = std(r_all_3_post_smkv_full);
s3_post_smkv2 = std(r_all_3_post_smkv);


% data normality
figure('DefaultAxesFontSize',12,'Position', [1 1 800 300])
subplot(1,2,1)
pd = fitdist(r_all_3_prior_smkv_full,distname); qqplot(r_all_3_prior_smkv_full, pd)
title('')
ylabel(sprintf('Quantiles of residuals prior Symkevi (%d)',length(r_all_3_prior_smkv_full) ))
grid('on')
subplot(1,2,2)
pd = fitdist(r_all_3_post_smkv_full,distname); qqplot(r_all_3_post_smkv_full, pd)
title('')
ylabel(sprintf('Quantiles of residuals during Symkevi (%d)',length(r_all_3_post_smkv_full) ))
grid('on')
saveas(gcf,fullfile(plotfolder,'fevModelBasedAnalysis_qqplot_none_smkv.png'))

[F_h_rejected_3, F_p_3, F_ci_3, F_stats_3] = vartest2(r_all_3_prior_smkv_full,r_all_3_post_smkv_full,'Tail','right','Alpha',alpha);

% Levene test
[L_h_rejected_3, L_p_3, L_stat_3, L_critical_3] = levenetest(r_all_3_prior_smkv_full,r_all_3_post_smkv_full,alpha);

%% function

function [residuals, curve] = applyMovingMean(x, y, w, t)
% t is inclusive
   curve = NaN(length(x),1);
   for i = 1:length(x)
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

function out = initPatientStatTable(n_rows);
    out = table('Size',[n_rows, 8],...
            'VariableTypes',["uint8","doublenan","doublenan","doublenan","doublenan","doublenan","doublenan","doublenan"],...
            'VariableNames',["p_processed_patients", "all_std", "all_n_residuals", "patientMaxVal", "std_residuals_prior_tripleT", "n_residuals_prior_tripleT", "std_residuals_post_tripleT","n_residuals_post_tripleT"]);
end

function [h, p, W, F] = levenetest(x1,x2,alpha)
    % upper-tail Levene test for homoscedasticity
    
    trimfactor = 20; % trimming factor, in matlab this corresponds to removing top 10% and bottom 10%
    k = 2;
    N1 = length(x1);
    N2 = length(x2);
    Z1 = abs(x1 - trimmean(x1,trimfactor)); % use median for Brown-Forsythe test
    Z2 = abs(x2 - trimmean(x2,trimfactor)); % use median for Brown-Forsythe test 
    N = N1 + N2;
    Z = cat(1,Z1,Z2);
    
    % (mean of the group residuals' absolute value - mean of all the
    % residuals' absolute value)^2, for each group
    num =  N1*(mean(Z1) - mean(Z))^2 + N2*(mean(Z2) - mean(Z))^2;
    
    % sum( (group residuals' absolute value - mean of the group residuals' 
    % absolute value)^2 ), for each group
    denom = sum( (Z1 - mean(Z1)).^2 ) + sum( (Z2 - mean(Z2)).^2 );
    
    % Levene ANOVA statistic
    W = (N-k)/(k-1) * num/denom;

    % critical threshold
    F = finv(1-alpha,k-1,N-k);
    
    % rejection
    h = W >= F;
    
    if h == 1 % find and approximation of the p-value
        while W >= finv(1-alpha,k-1,N-k)
            alpha = alpha/1.001;
        end
        p = alpha*1.001;
    else
         while W < finv(1-alpha,k-1,N-k)
            alpha = alpha*1.001;
         end
        p = alpha/1.001;
    end
end