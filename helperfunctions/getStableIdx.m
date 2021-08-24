% filters indexes corresponding to patient's stable days
% 
% Input:
% ------
% patient             -
% FEVdata             FEV data, uses serial date only
% treatments_table    ivandmeasures
% n_prior_t           #days prior to treatment start
% n_post_t            #days after modulator end
% modulators_table    brDrugTherapy
% npost_m             #days after modulator start
% filter              1 applies t and m, 2 applies t, 3 applies m

% Output: for selected patient
% -------
% idx2keep          stable days period - over whole FEVdata vector
% days_t            treatments' influence time window
% days_m            modulators' influence time window

function [idx2keep, days_t, days_m] = getStableIdx(patient, FEVdata, ...
    treatments_table, nprior_t, npost_t, modulators_table, npost_m, filter)

    % concatenate days to remove for treatments
    days_t=[];
    treatments_table = treatments_table(treatments_table.SmartCareID == patient,:); % faster
    for i = 1:size(treatments_table,1)
        days_t = [days_t (treatments_table.IVDateNum(i) - nprior_t) ...
            : (treatments_table.IVStopDateNum(i) + npost_t)];
    end
    days_t = unique(days_t); % remove duplicates
    
    % idem for modulators
    days_m=[];
    modulators_table = modulators_table(modulators_table.ID == patient, :);
    for i = 1:size(modulators_table,1)
        days_m = [days_m modulators_table.DateNum(i): ...
            modulators_table.DateNum(i) + npost_m];
    end
    days_m = unique(days_m); % remove duplicates
    
    % create mask with corresponding idx
    switch filter
        case 1 % treatments and modulators
            idx2keep = FEVdata(:,1) == patient ...
                & not(ismember(FEVdata(:,2), days_t)) ...
                & not(ismember(FEVdata(:,2), days_m));
        case 2 % only treatments
            idx2keep = FEVdata(:,1) == patient ...
                & not(ismember(FEVdata(:,2), days_t));
        case 3 % only modulators
            idx2keep = FEVdata(:,1) == patient ...
                & not(ismember(FEVdata(:,2), days_m));
    end
end