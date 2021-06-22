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

%% high level 

fprintf('amIntreNormcube contains %i values, %i NaN, %i zeros\n',numel(amIntrNormcube), sum(sum(sum(isnan(amIntrNormcube)))), sum(sum(sum(amIntrNormcube == 0))))

%% 

% TODO % issue with code from Damian, last day of align_wind is not normalised

% derive high level parameters
align_wind = size(amIntrNormcube,2)-1;% TODO % remove - 1

% initialise outputs
z = zeros(align_wind,length(measures.Index));
s = z;
n = z;

% do for one measure
m = 1; % calories

%% 

for m = measures.Index
    % get z, s for each point of the latent profile
    for d = 1:align_wind
        % over all interventions
         z(d,m) = mean(amIntrNormcube(:,d,m),'omitnan');
         s(d,m) = std(amIntrNormcube(:,d,m),'omitnan');
         n(d,m) = sum(~isnan(amIntrNormcube(:,d,m)))/size(amIntrNormcube,1)*100;
    end
end

%% plot

figure('DefaultAxesFontSize',12,'Position', [1 1 2000 1000])

t = tiledlayout(4,4);
meas = 1:17; meas(5)=[];

for m = meas %[1, 2, 6, 8, 12, 14, 16, 17]
    ax = nexttile;
    plotresults(ax,0:align_wind-1,z(:,m),s(:,m),n(:,m),measures.DisplayName{m},amIntrNormcube);
end

saveas(gcf,fullfile(plotfolder,'UnalignedModel_munorm_before_ape.png'))

function plotresults(ax,x,z,s,n,m_name,data)

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