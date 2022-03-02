function createMeasuresHeatmap(physdata, offset, cdPatient, study)

% createMeasuresHeatmap - creates the Patient/Measures heatmap

fprintf('Creating Heatmap of Measures\n');
fprintf('----------------------------\n');
tic

temp = hsv(64);
brightness = 0.9;


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
    colors(2,:)  = temp(5,:);
    colors(3,:)  = temp(6,:);
    colors(4,:)  = temp(7,:);
    colors(5,:)  = temp(8,:);
    colors(6,:)  = temp(9,:);
    colors(7,:)  = temp(10,:);
    colors(8,:)  = temp(11,:);
    colors(9,:)  = temp(12,:);
    colors(10,:)  = temp(13,:);
    colors(11,:)  = temp(14,:);
    colors(12,:)  = temp(15,:);
    colors(13,:)  = temp(16,:);
    colors(14,:)  = temp(17,:);
    colors(15,:)  = temp(18,:);
    colors(16,:)  = temp(19,:);
    colors(17,:)  = temp(20,:);
    colors(18,:)  = temp(21,:);
    colors(19,:)  = [1 0 1];
    nmeasures = 18;
else
    fprintf('**** Unknown Study ****');
    return;
end

colors(1:end - 1, :) = colors(1:end - 1,:) .* brightness;

% create a table of counts of measures by patient/day (@max function here
% is irrelevant as we just want the group counts
pdcountmtable = varfun(@max, physdata(:, {'SmartCareID','DateNum'}), 'GroupingVariables', {'SmartCareID', 'DateNum'});

% create the heatmap
title = sprintf('%s-Heatmap of Measures (%i patients, start date: %s, end date: %s)', study, ...
    length(unique(pdcountmtable.SmartCareID)), datestr(offset), datestr(offset + max(pdcountmtable.DateNum)));

[f, p] = createFigureAndPanel(title, 'portrait', 'a4');
h = heatmap(p, pdcountmtable, 'DateNum', 'SmartCareID', 'Colormap', colors, 'MissingDataColor', 'black', ...
    'ColorVariable','GroupCount','ColorMethod','max', 'MissingDataLabel', 'No data');
h.Title = ' ';
h.XLabel = 'Days';
h.YLabel = 'Patients';
%h.YLimits = {dispmin,dispmax};
h.ColorLimits = [0.5 size(colors,1)+0.5];
h.CellLabelColor = 'none';
h.GridVisible = 'off';

%[C,x] = sortx(h);

% save results
basedir = setBaseDir();
filename = sprintf('%s-HeatmapAllPatients', study);
subfolder = sprintf('Plots/%s', study);
if ~exist(strcat(basedir, subfolder), 'dir')
    mkdir(strcat(basedir, subfolder));
end
savePlotInDir(f, filename, subfolder);
close(f);

toc
fprintf('\n');

end
