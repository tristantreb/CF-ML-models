% plot the measures profile for each intervention used in the model
%
% 1) bar plot of interventions' duration
% 
% 
% Input:
% ------
% clinical and measurements data
% BRalignmentmodelinputs_recovery_gap*.mat
% BRmuNorm.mat          mu normalisation for mean window
%
% Output:
% -------
% one plot per patient

init;

modelinputfile = 'BRalignmentmodelinputs_recovery_gap10.mat';
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
% clean modulators tables
brDrugTherapy.DrugTherapyType = cleanDrugTherapyNamings(brDrugTherapy.DrugTherapyType);
% adds columns with serial date num
brDrugTherapy.DateNum = datenum(brDrugTherapy.DrugTherapyStartDate) - broffset;

%% meta information about interventions

figure('DefaultAxesFontSize',12,'Position', [1 1 1500 600])
subplot(2,1,1)
barInterventions(amInterventions.IVDateNum,amInterventions.IVStopDateNum,...
    sprintf('The %i interventions grouped by duration (post data completeness filter)',size(amInterventions,1)))

% bonus: same plot for all interventions
%load(fullfile(basedir, subfolder, 'BRivandmeasures_recovery_gap10.mat')); % note BRivandmeasures_gap*.mat also works

subplot(2,1,2)
barInterventions(ivandmeasurestable.IVDateNum,ivandmeasurestable.IVStopDateNum,...
    sprintf('The %i interventions grouped by duration (all interventions included)',size(ivandmeasurestable,1)))

%% plot
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
        
    for m = [14, 17, 6, 2, 3, 12, 16, 8]
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
    legend('Values',[amInterventions.Route{i} ' treatment'],'Meanwindow','Normmean')
    % write title
    sgtitle(sprintf('Intervention %i, patient %i', i, id))
    saveas(gcf,fullfile(plotfolder,sprintf('Intervention%i_ID%i.png', i, id)))
    close all;
end

%% functions
function barInterventions(start, stop, plottitle)
% bar plot with interventions grouped by duration
% Input: start and stop define the span of each intervention
% Output: bar plot

Intr = table(stop - start,'VariableNames',{'Duration'});
Intr = groupcounts(Intr, 'Duration');
bar(Intr.Duration, Intr.GroupCount)
grid('on')
title(plottitle)
xlabel('Intervention duration (days)')
ylabel('Frequency')
xticks(min(Intr.Duration):max(Intr.Duration))
yticks(0:10:max(Intr.GroupCount)+10)
end

function out = getRouteColor(route)
% assign the color for the corresponding route
switch route
    case 'Oral'
        out = [1     0.85  0   ];
    case {'IV', 'IVPBO'}
        out = [1     0     0   ];
end
end