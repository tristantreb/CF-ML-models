init;

% change during model result selection
nl = 1;
rs = 4;
prefix = 'tb_failnolist';

% changed during parameter optimisation
dw = 20;
mm = 34;
vm = 0.4;
ou = 15;

% load all files that are similar except their random seed
dircontent = dir(fullfile(basedir, subfolder, sprintf('%s*BRvEMMC_gp10_lm1_sig4_mu4_ca2_sm2_rm4_in1_im1_cm2_mm%i_od0_ou%i_dw%i_nl%i_rs%i_ds1_ct5_sc22-V_vs1_vm%.1f*.mat',prefix,mm, ou, dw, nl, rs, vm)));
ModelResultsFiles = cell(size(dircontent,1),1);
for a = 1:size(ModelResultsFiles,1)
    ModelResultsFiles{a} = dircontent(a).name;
end

%% load file data 
nfiles = size(ModelResultsFiles,1);
fprintf('Loading results\n');
fprintf('---------------\n');

for i = 1%:nfiles

    fprintf('Model number %3i - Filename: %s\n', i, ModelResultsFiles{i});
    load(fullfile(basedir, subfolder,ModelResultsFiles{i}),'meancurvemean', 'measures','measuresmask','amInterventions', 'offset','ninterventions'); 
end
initamInterventions = amInterventions;
initial_mean = meancurvemean;

%% set recovery start for all following plots

day_recovery_start = 10;
fontsize = 15;

if nl == 2
%     % assign C1 to largest group
%     idxC1 = amInterventions.LatentCurve == 1;
%     if sum(idxC1) < sum(not(idxC1))
%         amInterventions.LatentCurve(idxC1) = 2;
%         amInterventions.LatentCurve(not(idxC1)) = 1;
%     end
    
    idx=1;
    meancurvemean = initial_mean(idx,:,:);
    amInterventions = initamInterventions(initamInterventions.LatentCurve == idx,:);
    
end

%% response to treatment plot

figure('DefaultAxesFontSize',fontsize,'Position', [1 1 800 400])
tao = computeResponseToTreatment(amInterventions, day_recovery_start);
histogram(tao)
xlabel('\tau (days)');
ylabel(sprintf('Proportion of interventions (%%)', size(amInterventions,1)));
title(sprintf('mean = %.1f days (std = %.1f), median = %.0f days', mean(tao),std(tao),median(tao)));
saveas(gcf,fullfile(plotfolder,string([ModelResultsFiles{1}(1:end-4) '.png'])))


%% main plot (typical profile + bar graph with treatments start/end

% create figure

figure('DefaultAxesFontSize',fontsize,'Position', [1 1 800 800])

% define x axis
x = -day_recovery_start : (size(meancurvemean,2) - day_recovery_start -1 );


% plot typical profile
ax1 = subplot(4,1,[1 2 3]);
hold on

% define parameters
mm = 22; % 34 don't plot FEF2575 because very similar to FEV1 (helped convergence)
measures = getMeasuresMask(mm, measures);  measures.Mask(16) = false; % temperature not meaningful
% set method to align measures
align = "min";
align = "consensus_end";
align = "consensus_start";
align = "no";
% colors for measure - same as APE study
mcolors = [[0 0.4470 0.7410] ... %blue
    %;[0.9290 0.6940 0.1250]...%orange
    ;[0.8500 0.3250 0.0980]...%red
    ;[0.4660 0.6740 0.1880] ...%g
    ;[0.8 0.8 0] ... %yellow
    ;[0.4940 0.1840 0.5560] ...%purple
    %;[0.6350 0.0780 0.1840] ... %bordeau
    ;[0.3010 0.7450 0.9330]]; % light blue
color = 1;

% plot
idx = logical(measures.Mask); %idx(16)=false; 
for m = measures.Index(idx)'
    smoothedmeancurvemean = smartSmooth(meancurvemean(1,:,m), 3);
    y = adjustMeanCurve(smoothedmeancurvemean,align,day_recovery_start)';
%     if ismember(measures.DisplayName(m),getInvertedMeasures(study))
%         y = -y;
%     end
    if ismember(string(measures.DisplayName(m)),["Wellness", "FEV1"])
        plot(x,y,'LineStyle','-', 'Color', mcolors(color,:),'LineWidth',6,'Marker','.','MarkerSize',3);
    else
        plot(x,y,'LineStyle','-', 'Color', mcolors(color,:),'LineWidth',3,'Marker','.','MarkerSize',3);
    end
    yl = ylim;
    color = color+1;
end
%fill([-day_recovery_start 35-day_recovery_start 35-day_recovery_start -day_recovery_start], ...
%            [0 -0.27 -0.27 0], 'k', 'FaceAlpha', '0.1', 'EdgeColor', 'none');
hold off
grid('on');
ylabel('Signal change from stable baseline');
xlim([min(x),max(x)]);
if mod(min(x),2) == 0
    xticks(min(x):2:max(x));
else
    xticks(min(x)+1:2:max(x));
end

% annotations
xline(0,'--k','LineWidth',1);
x_annotation =  0.02 + 0.08 + (0.903-0.08)/34*day_recovery_start; % linear interpolation between max and min location of the graph window
annotation('textbox',[.9 .803 .1 .05],'String','stable baseline','FontSize',12,'FontWeight','bold','EdgeColor','none'); % faster by hand
yline(0,'--k','LineWidth',1);
annotation('textbox',[max(x_annotation,0.13) .85 .1 .1],'String','recovery start','FontSize',12,'FontWeight','bold','EdgeColor','none');
%legend

mask = ismember(measures.DisplayName,getInvertedMeasures(study));
measures.DisplayNameLegend(not(mask)) = measures.DisplayName(not(mask));
for i = 1:size(measures,1)
    if mask(i) == 1
        measures.DisplayNameLegend(i) = {measures.DisplayName{i} + " inverted"};
    end
    if i == 13
        measures.DisplayNameLegend(i) = {'O2 saturation'};
    elseif i == 11
        measures.DisplayNameLegend(i) = {'Minutes asleep'};
    end
end

LegendsStrings = measures.DisplayNameLegend(idx);
legend(LegendsStrings,'Interpreter','none','location','southeast');

%day_recovery_start

% plot bar graphs
ax2 = subplot(4,1,4); 

% treatment start
Date = amInterventions.Offset - day_recovery_start;
TreatmentStart = table(Date);
TreatmentStart = groupcounts(TreatmentStart,'Date');
bar(ax2, TreatmentStart.Date, TreatmentStart.GroupCount/sum(TreatmentStart.GroupCount)*100,'FaceColor',[0 0.2 0.12],'FaceAlpha',0.8,'EdgeColor','k','LineWidth',0.1)

% treatment end
hold on
Date = amInterventions.IVStopDateNum - amInterventions.IVDateNum + amInterventions.Offset - day_recovery_start;
TreatmentStop = table(Date);
TreatmentStop = groupcounts(TreatmentStop,'Date');
idx_in = TreatmentStop.Date < max(x);
n_interventions_stopped_after_xmax = sum(not(idx_in));
TreatmentStop = TreatmentStop(idx_in ,:);
bar(ax2, TreatmentStop.Date, TreatmentStop.GroupCount/sum(TreatmentStop.GroupCount)*100,'FaceColor',[.55 0 0],'FaceAlpha',0.8,'EdgeColor','k','LineWidth',0.1)
text = sprintf('(%.0f%% treatments end after day %i)',n_interventions_stopped_after_xmax/sum(TreatmentStop.GroupCount)*100, max(x));
annotation('textbox',[0.72 0.19 .1 .1],'String',text,'FontSize',10,'EdgeColor','none');

% common definitions for the two bar graphs
hold off
ylim([0 max(10,max(max(TreatmentStop.GroupCount),max(TreatmentStart.GroupCount)))])
if mod(min(x),2) == 0
    xticks(min(x):2:max(x));
else
    xticks(min(x)+1:2:max(x));
end
legend('Start','End');
ylabel('Treatments (%)');
grid('on')
set(gca,'Fontsize',fontsize)

% common x axis for typical profile and bar graphs
xlabel('Number of days from recovery start'); 

linkaxes([ax1,ax2],'x');

saveas(gcf,fullfile(plotfolder,sprintf('%sprofile_align_%s.png',prefix,align)))
%close all;

%% plot reference profile for one measure

day_recovery_start=7;
x = -day_recovery_start : (size(meancurvemean,2) - day_recovery_start -1 );
mm = 22; % 34 don't plot FEF2575 because very similar to FEV1 (helped convergence)
measures = getMeasuresMask(mm, measures);  measures.Mask(16) = false; % temperature not meaningful
% set method to align measures
align = "no";

% plot
idx = logical(measures.Mask); %idx(16)=false; 
for m = 6% measures.Index(idx)'
    smoothedmeancurvemean = smartSmooth(meancurvemean(1,:,m), 3);
    y = adjustMeanCurve(smoothedmeancurvemean,align,day_recovery_start)';
    plot(x,y,'LineStyle','-', 'Color', 'k','LineWidth',3,'Marker','.','MarkerSize',3);
    yl = ylim;
end
%fill([-day_recovery_start 35-day_recovery_start 35-day_recovery_start -day_recovery_start], ...
%            [0 -0.27 -0.27 0], 'k', 'FaceAlpha', '0.1', 'EdgeColor', 'none');
hold off
grid('on');
ylabel('Signal change from stable baseline');
xlim([min(x),max(x)]);
if mod(min(x),2) == 0
    xticks(min(x):2:max(x));
else
    xticks(min(x)+1:2:max(x));
end

% annotations
xline(0,'--k','LineWidth',1);
x_annotation =  0.03 + 0.08 + (0.903-0.08)/34*day_recovery_start; % linear interpolation between max and min location of the graph window
annotation('textbox',[.9 .843 .1 .05],'String','stable baseline','FontSize',12,'FontWeight','bold','EdgeColor','none'); % faster by hand
yline(0,'--k','LineWidth',1);
annotation('textbox',[max(x_annotation,0.13) .85 .1 .1],'String','recovery start','FontSize',12,'FontWeight','bold','EdgeColor','none');
%legend
xlabel('Number of days from recovery start');
set(gca,'Fontsize',fontsize)



%% functions
function curve = adjustMeanCurve(meancurvemean,align,day_recovery_start)
day_max_recovery = 25;
switch align
    case "no"
         curve = meancurvemean(:);
    case "consensus_end"
        curve = meancurvemean(:) - meancurvemean(day_max_recovery);
    case "consensus_start"
        curve = meancurvemean(:) - meancurvemean(day_recovery_start);
    case "min"
        curve = meancurvemean(:) - min(meancurvemean(:));
end

end

function tao = computeResponseToTreatment(amInterventions, day_start) 
ninterventions = size(amInterventions,1);
tao = nan(ninterventions,1);
for i = 1:ninterventions
    if amInterventions.Offset(i) < day_start
        tao(i) = day_start - amInterventions.Offset(i);
    else
        tao(i) = 0;
    end 
end
end

function curveout = smartSmooth(curvein,smoothparam)
% start smoothing when there is data
% nan values are kept to nan in curveout
idx = ~isnan(curvein);

curveout(idx) = smooth(curvein(idx),smoothparam);
curveout(not(idx)) = curvein(not(idx));
end