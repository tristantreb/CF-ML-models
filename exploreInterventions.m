% explores the list of interventions
%
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
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brDrugTherapy');

%% explore drug therapies
getDrugTherapyInfo(brDrugTherapy);

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
    
    measurestoplot = ["FEV1", "O2Saturation",...
            "FEF2575", "Weight",...
            "PulseRate","Wellness",...
            "RestingHR","Cough"];
    
    for m = mapMeasuresToIndex(measurestoplot,measures)
        
        nexttile;
        
        % get raw data
        data = getMeasureTable(brphysdata,measures.Name{m},measures.Column{m});
        % get data specific to current intervention
        data = data( ismember( data.SmartCareID, id ) & ...
            ismember( data.DateNum, range ), : );

        % plot
        y = eval(sprintf('data.%s',measures.Column{m}));
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
        %plot(data.DateNum-start, fillmissing(y,'movmean',10),'ro');

        yl = ylim;
        fill([0 stop-start stop-start 0], ...
            [yl(1) yl(1) yl(2) yl(2)], getRouteColor(amInterventions.Route{i}), 'FaceAlpha', '0.1', 'EdgeColor', 'none');
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
    legend('Values','Smoothed curve', [amInterventions.Route{i} ' treatment'],'Meanwindow','Normmean','Location','southwest')
    % write title
    sgtitle(sprintf('Intervention %i (from %s to %s), patient %i, smooth %i', i, datestr(broffset-1+start), datestr(broffset-1+stop), id, smoothing_factor))
    saveas(gcf,fullfile(plotfolder,sprintf('Intervention%i_ID%i.png', i, id)))
    close all;
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