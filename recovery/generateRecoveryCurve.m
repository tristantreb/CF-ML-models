% finds the typical profile of the recovery with zero offset
%
% Input:
% ------
% *intrnormdatacube_recovery.mat
%
% Output:
% -------
%

init;
filename = sprintf('%sintrnormdatacube_recovery.mat', study);
fprintf('Loading alignment model Inputs data %s\n', filename);
load(fullfile(basedir, subfolder, filename));

% data window analysis

% size, #Nan, #0
fprintf('amIntreNormcube contains %i values, %i NaN, %i zeros\n',numel(amIntrNormcube), sum(sum(sum(isnan(amIntrNormcube)))), sum(sum(sum(amIntrNormcube == 0))))

% TODO % issue with code, last day of align_wind is not normalised

% align_wind
align_wind = size(amIntrNormcube,2)-1;% TODO % remove - 1
ninterventions = size(amIntrNormcube,1);

%% number of points
figure('DefaultAxesFontSize',12,'Position', [1 1 2000 1000])
tiledlayout(4,4);
meas = 1:17; meas(5)=[];
for m = meas
    nexttile;
    barHistogram(sum(~isnan(amIntrNormcube(:,1:align_wind,m)),2)',...
        measures.DisplayName{m},...
        'Days with records');
end
title=sprintf('#recordings over the %i days data window (%i selected interventions)', align_wind, size(amIntrNormcube,1));
sgtitle(title);
saveas(gcf,fullfile(plotfolder,[title '.png']));
close all;

%% span of recording

figure('DefaultAxesFontSize',12,'Position', [1 1 2000 1000])
tiledlayout(4,4);
for m = meas
    span=nan(ninterventions,1);
    nexttile;
    %find position of first and last non nan element
    for i=1:ninterventions
        logic=~isnan(amIntrNormcube(i,1:align_wind,m));
        idx=find(logic==1);
        if length(idx) >= 1
            span(i) = idx(end)-idx(1)+1;
        else
            span(i)=0;
        end
    end
    barHistogram(span,...
        measures.DisplayName{m},...
        'Days with records')
end
title=sprintf('Span of recordings per measure (%i days data window, %i selected interventions)', align_wind, size(amIntrNormcube,1));
sgtitle(title);
saveas(gcf,fullfile(plotfolder,[title '.png']));
close all;

%%

% initialise outputs
z = zeros(align_wind,length(measures.Index));
s = z;
n = z;

% calculate z, s, n
for m = measures.Index
    % get z, s for each point of the latent profile
    for d = 1:align_wind
        % over all interventions
         z(d,m) = mean(amIntrNormcube(:,d,m),'omitnan');
         s(d,m) = std(amIntrNormcube(:,d,m),'omitnan');
         n(d,m) = sum(~isnan(amIntrNormcube(:,d,m)))/size(amIntrNormcube,1)*100;
    end
end

%% plot #measures latent curves

figure('DefaultAxesFontSize',12,'Position', [1 1 2000 1000])

tiledlayout(4,4);
meas = 1:17; meas(5)=[];

for m = meas
    ax = nexttile;
    plotresults(ax,0:align_wind-1,smooth(z(:,m),3),smooth(s(:,m),3),n(:,m),measures.DisplayName{m},amIntrNormcube);
end

saveas(gcf,fullfile(plotfolder,'UnalignedModel_munorm_before_ape_40d.png'))

%% violently superpose all interventions curves for selected measure 
% not conclusive

figure('DefaultAxesFontSize',12,'Position', [1 1 2000 1000])
tiledlayout(4,2);
meas = 1:17; meas(5)=[];
for m = [14, 17, 6, 2, 3, 12, 16, 8] %meas
    nexttile;
    hold on
    for i = 1:size(amIntrNormcube,1)
        plot(1:align_wind,smooth(amIntrNormcube(i,1:end-1,m),3))
    end
    ylabel(sprintf('%s Normalised',measures.DisplayName{m}))
    hold off
end
sgtitle(sprintf('Violently superpose all %i interventions curves',size(amIntrNormcube,1)))
saveas(gcf,fullfile(plotfolder,'All_timeseries_superposition.png'))

%% curve clustering

plot(smooth(fillmissing(amIntrNormcube(35,1:end-1,6),'linear'),5))


%% functions

function plotresults(ax,x,z,s,n,m_name,data)
% plot latent curve for one measure
yyaxis(ax,'left')
plot(x,z,'.b',x,z-s,'.b',x,z+s,'.b')
hold on
plot(x,z,'b',x,z-s,'--b',x,z+s,'--b')
hold off
ax.YColor = 'b';
ylabel(ax,sprintf('%s Normalised',m_name))

yyaxis(ax,'right')
bar(x,n,'k','FaceAlpha',0.1,'EdgeAlpha',0.0)
ax.YColor = 'k';
ylim(ax,[0,size(data,1)*4])
yticks(ax,[0 50 100])
ylabel(ax,'Data count (%)')

% ax.Color, background color
% ax.YColor, y axis color
% ax1.Box = 'off', remove ticks on the other side

xlabel('Days post treatment')

end