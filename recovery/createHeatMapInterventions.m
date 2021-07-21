init;

% load data brphysdata, broffset, brPatient, study
[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[physdata, offset, ~] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);
load(fullfile(basedir, subfolder, 'breatheclinicaldata.mat'),'brPatient');
cdPatient=brPatient; clear brPatient;

modelinputfile = 'BRalignmentmodelinputs_recovery_gap10_datawind_20.mat';
% load amInterventions, amDatacube, measures and count info
load(fullfile(basedir, subfolder, modelinputfile));

temp = hsv(64);
brightness = 0.9;

%% 
if ismember(study, {'SC', 'TM'})
    colors(1,:)  = temp(4,:);
    colors(2,:)  = temp(6,:);
    colors(3,:)  = temp(8,:);
    colors(4,:)  = temp(10,:);
    colors(5,:)  = temp(12,:);
    colors(6,:)  = temp(14,:);
    colors(7,:)  = temp(16,:);
    colors(8,:)  = temp(18,:);
    colors(9,:)  = temp(20,:);
    colors(10,:)  = [1 0 1];
    nmeasures = 9;
elseif ismember(study, {'CL'})
    colors(1,:)  = temp(4,:);
    colors(2,:)  = temp(6,:);
    colors(3,:)  = temp(7,:);
    colors(4,:)  = temp(8,:);
    colors(5,:)  = temp(9,:);
    colors(6,:)  = temp(10,:);
    colors(7,:)  = temp(11,:);
    colors(8,:)  = temp(12,:);
    colors(9,:)  = temp(13,:);
    colors(10,:)  = temp(14,:);
    colors(11,:)  = temp(15,:);
    colors(12,:)  = temp(16,:);
    colors(13,:)  = temp(17,:);
    colors(14,:)  = temp(18,:);
    colors(15,:)  = temp(20,:);
    colors(16,:)  = [1 0 1];
    nmeasures = 15;
elseif ismember(study, {'BR'})
    colors(1,:)  = temp(4,:);
    colors(2,:)  = temp(6,:);
    colors(3,:)  = temp(7,:);
    colors(4,:)  = temp(8,:);
    colors(5,:)  = temp(9,:);
    colors(6,:)  = temp(10,:);
    colors(7,:)  = temp(11,:);
    colors(8,:)  = temp(12,:);
    colors(9,:)  = temp(13,:);
    colors(10,:)  = temp(14,:);
    colors(11,:)  = temp(15,:);
    colors(12,:)  = temp(16,:);
    colors(13,:)  = temp(17,:);
    colors(14,:)  = temp(18,:);
    colors(15,:)  = temp(20,:);
    colors(16,:)  = temp(21,:);
    colors(17,:)  = temp(22,:);
    colors(18,:)  = [1 0 1]; % purple for treatment end
    colors(19,:)  = temp(32,:); % blue for treatment start
    nmeasures = 19;
else
    fprintf('**** Unknown Study ****');
    return;
end

colors(1:end , :) = colors(1:end,:) .* brightness;

% % get the date scaling offset for each patient
% patientoffsets = getPatientOffsets(physdata);

%% 
days_prior = 35;
days_post = 39;
intrwindow = days_prior+days_post+1;
intrrangetable = table('Size',[ninterventions*intrwindow 6],...
    'VariableTypes',{'double','double','double','string','string','double'},...
    'VariableNames',{'ID','Date','IntrDate','Label','Route','Index'});
% nan(ninterventions*intrwindow,4);

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

% create a table of counts of measures by patient/day (@max function here
% is irrelevant as we just want the group counts
pdcountmtable = varfun(@max, physdata(:, {'SmartCareID','DateNum'}), 'GroupingVariables', {'SmartCareID', 'DateNum'});

intrrangetable = outerjoin(intrrangetable,pdcountmtable,'Type', 'Left','Leftkeys',{'ID' 'Date'},'RightKeys',{'SmartCareID' 'DateNum'});

% if stop date belongs to data range, change color
% ncqnnot vectorize as many dates can be spotted for one intervention
for i = 1:size(intrrangetable,1)
    % ID and date equality
    if sum(ismember(amInterventions.SmartCareID,intrrangetable.ID(i)) ...
            & ismember(amInterventions.IVStopDateNum,intrrangetable.Date(i)))
        
        intrrangetable.GroupCount(i) = 18;
    elseif sum(ismember(amInterventions.SmartCareID,intrrangetable.ID(i)) ...
            & ismember(amInterventions.IVDateNum,intrrangetable.Date(i)))
        
        intrrangetable.GroupCount(i) = 19;
    end
end

% create the heatmap
title = sprintf('%s-Heatmap of Measures per Intervention', study);

f = figure('Name', title,'Position', [1 1 2000 1000]);
p = uipanel('Parent',f,'BorderType','none', 'BackgroundColor', 'white');
h = heatmap(p, intrrangetable, 'IntrDate', 'Label', 'Colormap', colors, 'MissingDataColor', 'black', ...
    'ColorVariable','GroupCount','ColorMethod','max', 'MissingDataLabel', 'No data');
h.Title = sprintf('Heatmap of measures for the %i interventions. Legend: cyan = start date, pink = stop date',ninterventions);
h.XLabel = 'Days';
h.YLabel = 'Intervention';
h.CellLabelColor = 'none';
h.GridVisible = 'off';

filename = sprintf('%s-HeatmapInterventions', study);
subfolder = sprintf('Plots/%s', study);
if ~exist(strcat(basedir, plotfolder), 'dir')
    mkdir(strcat(basedir, plotfolder));
end
set(gca,'FontSize',6);
saveas(gcf,fullfile(plotfolder,[filename '.png']))
close(f);