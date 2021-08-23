% explores the list of interventions
% 
% - plot patient enrolment time
% - plot number of interventions par patient
% - explore drug therapies (display two tables)
% - display bar plots of interventions durations
% - plot the measures profile for each intervention used in the model
% 
% 
% 
% Input:
% ------
% clinical and measurements data
% BRalignmentmodelinputs_recovery_gap*.mat  uses amInterventions ID, start,
% stop date, 
% BRmuNorm.mat                              mu normalisation for mean window
%
% Output:
% -------
% one plot per patient

init;

modelinputfile = 'BRalignmentmodelinputs_recovery_gap10_datawind20.mat';
munormfile = 'BRmuNorm.mat';

% load amInterventions, amDatacube, measures and count info
load(fullfile(basedir, subfolder, modelinputfile));
% load normmean
load(fullfile(basedir, subfolder, munormfile));

% load measures
[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[brphysdata, broffset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);

% load CFTR modulators therapy
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brDrugTherapy', 'brPatient');

load('/Users/tristan.trebaol/Documents/PDM/Project/MatlabSavedVariables/BRivandmeasures_recovery_gap10.mat')

%% explore drug therapies
getDrugTherapyInfo(brDrugTherapy, brPatient);

%% patient enrollment time

% patient enrolment time
newest_date = max(datenum(brPatient.PatClinDate));
idx = ~isnat(brPatient.WithdrawalDate);
brPatient.TimeEnrolled(idx) = datenum(brPatient.WithdrawalDate(idx)) - datenum(brPatient.StudyDate(idx));
brPatient.TimeEnrolled(not(idx)) = newest_date - datenum(brPatient.StudyDate(~idx));
brPatient.TimeEnrolled=brPatient.TimeEnrolled/365.25*12; % in months
figure('DefaultAxesFontSize',16,'Position', [1 1 600 300])
histogram(brPatient.TimeEnrolled,'BinWidth',3,'FaceColor','k','FaceAlpha',0.4,'EdgeColor','k','LineWidth',1)
xticks(0:3:30)
xlabel('Enrolment time (months)')
ylabel('Number of participants')
saveas(gcf,fullfile(plotfolder,sprintf('Patient_enrolment_time_from%s_to%s.png',min(brPatient.PatClinDate),max(brPatient.PatClinDate))));
close all;


%% display #interventions per patient

IDs = brPatient.ID(brPatient.TimeEnrolled > 0);
%nintr = groupcounts(amInterventions,'SmartCareID');
nintr = groupcounts(ivandmeasurestable,'SmartCareID');
% patient enrolled since > 12 months
nintr_patient = nintr(ismember(nintr.SmartCareID,IDs),:);
fprintf('Proportion of patients with at least one intervention %.2f\n', size(nintr,1)/size(IDs,1));

nintr = groupcounts(nintr_patient,'GroupCount');
nintr(end+1,:) = table( 0, size(unique(brPatient.ID),1) - size(unique(ivandmeasurestable.SmartCareID),1), 0);

figure('DefaultAxesFontSize',16,'Position', [1 1 600 300])
b = bar(nintr.GroupCount,nintr.GroupCount_1,'FaceColor','k','FaceAlpha',0.4,'EdgeColor','k','LineWidth',1);
xlabel('Number of antibiotic treatments');
ylabel('Number of participants');
b.FaceColor = 'flat'; b.CData(1,:) = [1 1 1];
saveas(gcf,fullfile(plotfolder,sprintf('bar_intrperpatient_before_filter.png')));
close all;

%% display bar plots of interventions durations

figure('DefaultAxesFontSize',12,'Position', [1 1 1500 600])
subplot(2,1,1)
barHistogram(amInterventions.IVStopDateNum-amInterventions.IVDateNum,...
    sprintf('The %i interventions grouped by duration (post data completeness filter)',size(amInterventions,1)),...
    'Intervention duration (days)')

% bonus: same plot for all interventions
load(fullfile(basedir, subfolder, 'BRivandmeasures_recovery_gap10.mat')); % note BRivandmeasures_gap*.mat also works

subplot(2,1,2)
barHistogram(ivandmeasurestable.IVStopDateNum-ivandmeasurestable.IVDateNum,...
    sprintf('The %i interventions grouped by duration (all interventions included)',size(ivandmeasurestable,1)),...
    'Intervention duration (days)')

saveas(gcf,fullfile(plotfolder,sprintf('Bar graph of the interventions duration (%i and %i).png', size(amInterventions,1), size(ivandmeasurestable,1))));
close all;

%% plot the measures profile for each intervention used in the model

% parameters
days_prior = 35; % include mu normalisation window [-35, -25] days
days_post = 39; % treatment generally durate 2 weeks, includes day 0

% TODO % add which Drug Therapy has been started, and if within range, plot it

% note amInterventions date 0 is study start date, i.e. broffset (not patient start date)

for i = 1:ninterventions
    figure('DefaultAxesFontSize',12,'Position', [1 1 1500 600])
    t = tiledlayout(4,2);
    
    id = amInterventions.SmartCareID(i);
    start = amInterventions.IVDateNum(i);
    stop = amInterventions.IVStopDateNum(i);
    range = start - days_prior : start + days_post;
    
    % visualize all treatments within range before or after
    start_in_range = amInterventions.IVDateNum(amInterventions.SmartCareID == id);
    stop_in_range = amInterventions.IVStopDateNum(amInterventions.SmartCareID == id);
    route_in_range = string(amInterventions.Route(amInterventions.SmartCareID == id));
    start_in_range(not(ismember(start_in_range, range))) = nan;
    stop_in_range(not(ismember(stop_in_range, range))) = nan;

    idx = isnan(start_in_range) & not(isnan(stop_in_range));
    start_in_range(idx) = range(1);
    idx = isnan(stop_in_range) & not(isnan(start_in_range));
    stop_in_range(idx) = range(end);
    
    start_in_range(isnan(start_in_range))=[];
    route_in_range(isnan(stop_in_range)) = [];
    stop_in_range(isnan(stop_in_range))=[];

    
    measurestoplot = ["FEV1", "Wellness",...
            "FEF2575", "Cough",...
            "PulseRate","O2Saturation",...
            "Temperature","MinsAsleep"];
%     measurestoplot = ["FEV1", "O2Saturation",...
%             "PulseRate", "HasColdOrFlu",...
%             "MinsAsleep","Wellness",...
%             "MinsAwake","Temperature"];
    
    for m = mapMeasuresToIndex(measurestoplot,measures)
        
        nexttile;
        
        % get raw data
        data = getMeasureTable(brphysdata,measures.Name{m},measures.Column{m});
        % get data specific to current intervention
        data = data( ismember( data.ID, id ) & ...
            ismember( data.DateNum, range ), : );

        % plot
        y = eval(sprintf('data.%s',measures.Name{m}));
        plot(data.DateNum-start, y,...
            'Color', [0, 0.65, 1], ...
            'LineStyle', ':', ...
            'Marker', 'o', ...
            'LineWidth',1,...
            'MarkerSize',2,...
            'MarkerEdgeColor','b',...
            'MarkerFaceColor','g');
        hold on
        smoothing_factor=5;
        plot(data.DateNum-start, smooth(y,smoothing_factor),...
            'Color', [0, 0.65, 1], ...
            'LineStyle', '-',...
            'LineWidth',1);
        if length(y)>1 && min(y) ~= max(y)
            ylim([min(y), max(y)]);
        end
        yl = ylim;
        for j = 1:length(start_in_range)
            startdate = start_in_range(j) - start;
            stopdate = stop_in_range(j) - start;
            route = route_in_range(j);
            fill([startdate stopdate stopdate startdate], ...
                [yl(1) yl(1) yl(2) yl(2)], getRouteColor(route), 'FaceAlpha', '0.1', 'EdgeColor', 'none');
        end
        fill([-35 -25 -25 -35], ...
            [yl(1) yl(1) yl(2) yl(2)], 'k', 'FaceAlpha', '0.04', 'EdgeColor', 'none');

        yline(normmean(i,m))
        hold off

        xlim([-days_prior days_post+1])
        xlabel('Days from study start')
        ylabel(measures.DisplayName{m})

        % reverse values where necessary
        if ismember(measures.DisplayName{m},getInvertedMeasures(study))
            set(gca, 'YDir','reverse')
        end
    end
    if length(route_in_range) == 1
        legend('Values','Smoothed curve', [amInterventions.Route{i} ' treatment'],'Meanwindow','Normmean','Location','eastoutside')
    elseif length(route_in_range) == 2
        if start_in_range(1) == start
            legend('Values','Smoothed curve', [amInterventions.Route{i} ' treatment'], [amInterventions.Route{i+1} ' treatment'],'Meanwindow','Normmean','Location','eastoutside')
        else
            legend('Values','Smoothed curve', [amInterventions.Route{i-1} ' treatment'], [amInterventions.Route{i} ' treatment'],'Meanwindow','Normmean','Location','eastoutside')
        end
    else
        continue;
    end
    % write title
    sgtitle(sprintf('Intervention %i (from %s to %s), patient %i, smooth %i', i, datestr(broffset-1+start), datestr(broffset-1+stop), id, smoothing_factor))
    saveas(gcf,fullfile(plotfolder,sprintf('Intervention%i_ID%i.png', i, id)))
    close all;
end

%%
mm = 22; % 34
measures = getMeasuresMask(mm, measures);
measures = measures.Index(logical(measures.Mask));
nmeasures=length(measures);

align_wind = 20;
i = 1;

r = nan(nmeasures, nmeasures);

for x = 1:measures
    for y = 1:nmeasures
        X = amIntrDatacube(i,1:align_wind,x)
        r(x,y) = xcorr(amIntrDatacube(i,1:align_wind,x),amIntrDatacube(i,1:align_wind,y));
    end
end
%% functions

function out = getRouteColor(route)
% assign the color for the corresponding route
switch route
    case 'Oral'
        out = [1     0.85  0   ];
    case {'IV', 'IVPBO'}
        out = [1     0     0   ];
end
end

function m_idx = mapMeasuresToIndex(m, measures)
% subsitutes the measure name by its index thans to the "measures" array
%
% Note: vectorisation is agnostic of the ordering of m, it takes the 
% ordering of measures. Hence we made this function
m_idx = zeros(1,length(m));
for i = 1:length(m)
    m_idx(1,i) = measures{ismember(measures.DisplayName, m(i)),'Index'};
end
end