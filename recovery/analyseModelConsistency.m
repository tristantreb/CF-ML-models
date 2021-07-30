init;

% differentiating parameters
nl = 2;
dw = 20;
mm = 22;

% load all files that are similar except their random seed
dircontent = dir(fullfile(basedir, subfolder, sprintf('BRvEMMC_gp10_lm1_sig4_mu4_ca2_sm2_rm4_in1_im1_cm2_mm%i_mo1_dw%i_nl%i_rs*_ds1_ct5_sc22-V_vs0_vm0.0*.mat',mm,dw, nl)));
ModelResultsFiles = cell(size(dircontent,1),1);
for a = 1:size(ModelResultsFiles,1)
    ModelResultsFiles{a} = dircontent(a).name;
end


nfiles = size(ModelResultsFiles,1);
fprintf('Loading results\n');
fprintf('---------------\n');

for i = 1:nfiles
    % link #column with filename
    fprintf('Model number %3i - Filename: %s\n', i, ModelResultsFiles{i});
    load(fullfile(basedir, subfolder,ModelResultsFiles{i}),'amInterventions');
    % add latent cruve allocation
    curveallocation(:,i+1) = amInterventions.LatentCurve;
end
curveallocation(:,1) = amInterventions.IntrNbr;
ninterventions = size(curveallocation,1);

% -> started to add rs in top row
% load(fullfile(basedir, subfolder,ModelResultsFiles{1}),'amInterventions');
% ninterventions = size(amInterventions,1);
% curveallocation = nan(ninterventions+1,nfiles+1);
% curveallocation(2:end,1) = amInterventions.IntrNbr;
% 
% 
% for i = 1:nfiles
%     % link #column with filename
%     fprintf('Model number %3i - Filename: %s\n', i, ModelResultsFiles{i});
%     load(fullfile(basedir, subfolder,ModelResultsFiles{i}),'amInterventions');
%     % add latent cruve allocation
%     curveallocation(2:end,i+1) = amInterventions.LatentCurve;
%     temp = split(ModelResultsFiles{i},'_');
%     %curveallocation(1,i+1) = str2num(temp{16});
% end


%% harmonise the class numbers bw comparing the columns

% row: class number
% column: random seed
harm = zeros(nl,nfiles+1);
harm(:,1) = 1:nl;
for i = 1:nfiles
    harm(:,i+1) = groupcounts(curveallocation(:,i+1));
end
fprintf('Result of the curve allocation for %i interventions:\n',ninterventions);
disp(harm)

%% sorted view

% TODO % should remember rs number (columns) and curve index (rows)

harmsorted = sort(harm,1);
[temp, order] = sort(harmsorted(1,:));
disp(harmsorted(:,order));

%% consistency of the curve allocation for different random seeds
accuracy = zeros(ninterventions,nl);
for i = 1:ninterventions
    accuracy(i,:) = getCurveAllocation(curveallocation(i,2:end),nl);
end

%% function


function out = getCurveAllocation(row,nl)
% return number of latent curves
out = zeros(1,nl);
for i = 1:nl
    for j = row
        if i == j
            out(i) = out(i)+1;
        end
    end
    %out(i) = out(i)*100/length(row);
end

end