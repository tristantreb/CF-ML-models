init;

% load data brphysdata, broffset, brPatient, study
[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[physdata, offset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brPatient');
cdPatient=brPatient; clear brPatient;

modelinputfile = 'BRalignmentmodelinputs_recovery_gap10_datawind20.mat';
% load amInterventions, amDatacube, measures and count info
load(fullfile(basedir, subfolder, modelinputfile));

%% 
% choose measure mask
daysprior = 35;
dayspost = 39;

for measuresmask = 34
    plotHeatMap(measuresmask, daysprior, dayspost, amInterventions, physdata, measures, study);
    % save value
    filename = sprintf('%s-HeatmapInterventions-mm%i', study, measuresmask);
    subfolder = sprintf('Plots/%s', study);
    if ~exist(strcat(basedir, plotfolder), 'dir')
        mkdir(strcat(basedir, plotfolder));
    end
    saveas(gcf,fullfile(plotfolder,[filename '.png']))
    close all;
end


%% helper functions

function plotHeatMap(mm, days_prior, days_post, amInterventions, physdata, measures, study)

intrwindow = days_prior+days_post+1;
ninterventions = size(amInterventions,1);

intrrangetable = table('Size',[ninterventions*intrwindow 6],...
    'VariableTypes',{'double','double','double','string','string','double'},...
    'VariableNames',{'ID','Date','IntrDate','Label','Route','Index'});

for i = 1:ninterventions
    start = amInterventions.IVDateNum(i);
    stop = amInterventions.IVStopDateNum(i);
    range = start - days_prior : start + days_post;
    
    intrrangetable.ID(intrwindow*(i-1)+1:intrwindow*i)=amInterventions.SmartCareID(i);
    intrrangetable.Date(intrwindow*(i-1)+1:intrwindow*i)=range;
    intrrangetable.IntrDate(intrwindow*(i-1)+1:intrwindow*i)=-days_prior:days_post;
    intrrangetable.Index(intrwindow*(i-1)+1:intrwindow*i)=i;
    intrrangetable.Route(intrwindow*(i-1)+1:intrwindow*i)=amInterventions.Route{i};
end

% add label
idxlabel = strsplit(sprintf(" - Intr %i - ,",intrrangetable.Index),',')';
idlabel = strsplit(sprintf("ID %i,",intrrangetable.ID),',')';
intrrangetable.Label = intrrangetable.Route + idxlabel(1:end-1) + idlabel(1:end-1);

% fitler corresponding recordings
measures = getMeasuresMask(mm, measures);
physdatafiltered = physdata(ismember(physdata.RecordingType,measures.Name(logical(measures.Mask))),:);

% create a table of counts of measures by patient/day (@max function here
% is irrelevant as we just want the group counts)
pdcountmtable = varfun(@max, physdatafiltered(:, {'SmartCareID','DateNum'}), 'GroupingVariables', {'SmartCareID', 'DateNum'});
% add groupcount to intervention table
intrrangetable = outerjoin(intrrangetable,pdcountmtable,'Type', 'Left','Leftkeys',{'ID' 'Date'},'RightKeys',{'SmartCareID' 'DateNum'});

% define colors
nmeasures = max(intrrangetable.GroupCount);
colors = getColors(nmeasures);

% if stop date belongs to data range, change color
% ncqnnot vectorize as many dates can be spotted for one intervention
for i = 1:size(intrrangetable,1)
    % ID and date equality
    if sum(ismember(amInterventions.SmartCareID,intrrangetable.ID(i)) ...
            & ismember(amInterventions.IVStopDateNum,intrrangetable.Date(i)))
        
        intrrangetable.GroupCount(i) = nmeasures+1;
    elseif sum(ismember(amInterventions.SmartCareID,intrrangetable.ID(i)) ...
            & ismember(amInterventions.IVDateNum,intrrangetable.Date(i)))
        
        intrrangetable.GroupCount(i) = nmeasures+2;
    end
end

% create the heatmap
title = sprintf('%s-Heatmap of Measures per Intervention', study);

f = figure('Name', title,'Position', [1 1 2000 1000]);
p = uipanel('Parent',f,'BorderType','none', 'BackgroundColor', 'white');
h = heatmap(p, intrrangetable, 'IntrDate', 'Label', 'Colormap', colors, 'MissingDataColor', 'black', ...
    'ColorVariable','GroupCount','ColorMethod','max', 'MissingDataLabel', 'No data');
h.Title = [sprintf('Heatmap of #measures for the %i interventions, mm = %i\n Measures:  %s\n',ninterventions,mm, ...
    string(join(measures.DisplayName(logical(measures.Mask)),', '))),...
    'Legend: {\color[rgb]{0.9 0 0.9}Interventions start}, {\color[rgb]{0 0.9 0.8156}Interventions end}'];

h.XLabel = 'Days';
h.YLabel = 'Intervention';
h.CellLabelColor = 'none';
h.GridVisible = 'off';
h.ColorLimits = [1 size(colors,1)+1];
set(gca,'FontSize',10);
end

function colors = getColors(nmeasures)
% create color scale for the heatmap
temp = hsv(64);
brightness = 0.9;

colors = zeros(nmeasures+2, 3);

idx = ceil(linspace(4,23,nmeasures));

colors(1:nmeasures,:) = temp(idx,:);
colors(end,:)  = [1 0 1]; % pink for treatment end
colors(end-1,:)  = temp(32,:); % cyan for treatment start

colors = colors .* brightness;
end