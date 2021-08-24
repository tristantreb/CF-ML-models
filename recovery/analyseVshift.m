% used to find the optimal maximum vertical shift value
% 
% Input: alignment results mat file
% ------
% 
% Output: histogram of the vertical shift's distribution 
% -------

init;

% differentiating parameters
nl = 1;
dw = 20;
mm = 34;
%vm=0.7;

% load all files that are similar except their random seed
dircontent = dir(fullfile(basedir, subfolder, sprintf('*BRvEMMC_gp10_lm1_sig4_mu4_ca2_sm2_rm4_in1_im1_cm2_mm%i_od0_ou15_dw%i_nl%i_rs*_ds1_ct5_sc22-V_vs1_vm*.mat',mm,dw, nl)));
ModelResultsFiles = cell(size(dircontent,1),1);
for a = 1:size(ModelResultsFiles,1)
    ModelResultsFiles{a} = dircontent(a).name;
end

nfiles = size(ModelResultsFiles,1);
fprintf('Loading results\n');
fprintf('---------------\n');

figure('DefaultAxesFontSize',17,'Position', [1 1 1600 400])
t = tiledlayout(1,3);
for i = 1:nfiles
    % link #column with filename
    fprintf('Model number %3i - Filename: %s\n', i, ModelResultsFiles{i});
    load(fullfile(basedir, subfolder,ModelResultsFiles{i}),'vshift', 'amInterventions','ninterventions','vshiftmax'); 
    if ismember(vshiftmax,[0.2 0.4 0.7])
        plotverticalshift(vshift,amInterventions.Offset,ninterventions,vshiftmax);
    end
end

%%
function plotverticalshift(vshift,offsets,ninterventions,vshiftmax)

% plot vshift for best matching offset
for i = 1:ninterventions
    best(i) = vshift(1,i,6,offsets(i)+1);
end

nexttile;
histogram(best',2*vshiftmax*10+1,'BinWidth',0.1)
title(['\alpha' sprintf(' = %.1f',vshiftmax)]);
ylim([0,55])
xlabel('Vertical shift of most probable offset (after convergence)')
ylabel('#interventions')
end