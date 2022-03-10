function [Drugsbypatients, historytable] = getDrugTherapyInfo(brDrugTherapy, brPatient) 
% all information about CFTR modulator therapies is centralised here
% 
% - groupcount the different drug therapies over patients
% - define drugs mix: list of the drug therapies given to a patient, sorted in time
% - groupcount the different drugs mix over patients and details the concerned patients IDs
% 
% Input: needs ID, DrugTherapyType, DrugTherapyStopDate from brDrugTherapy
% ------
% 
% Output: 
% -------
% - displays two tables
% - table with patients and their associated drug mix

% sort by date and ID
brDrugTherapy = sortrows(brDrugTherapy,'DrugTherapyStartDate','ascend');
brDrugTherapy = sortrows(brDrugTherapy,'ID','ascend');

patients = unique(brDrugTherapy.ID);

% group drug therapies by patients
History = string(zeros(length(patients),1));
Current = string(zeros(length(patients),1));
for i = 1:length(patients)
    drugs_list = brDrugTherapy(brDrugTherapy.ID == patients(i),{'DrugTherapyStopDate', 'DrugTherapyType', 'DrugTherapyStartDate'});
    
    % merge therapies that are consecutiveley the same 
    % TODO % this should be handled at source (in redcap)
    for j = 1:size(drugs_list,1)-1
        if strcmp(drugs_list.DrugTherapyType(j), drugs_list.DrugTherapyType(j+1))
            drugs_list.DrugTherapyStartDate(j+1) = drugs_list.DrugTherapyStartDate(j);
            drugs_list(j,:)=[];
        end
    end
    
    History(i,1) = join(string(drugs_list.DrugTherapyType),', ');
    Current(i,1) = string(drugs_list.DrugTherapyType(end));
    if ~isnat(drugs_list.DrugTherapyStopDate(end))
        History(i,1) = append(History(i,1),", _Therapy Stopped_");
        Current(i,1) = "_Therapy Stopped_";
    end
end

% (out var) table with patients and their associated history
Drugsbypatients = table(patients,History,Current);
Drugsbypatients = outerjoin(Drugsbypatients, brPatient, 'LeftKeys', 'patients', 'RightKeys', 'ID', 'LeftVariables', {'patients', 'History', 'Current'}, 'RightVariables', 'ID');
idxtoreplace = isnan(Drugsbypatients.patients);
Drugsbypatients.History(idxtoreplace) = "_No Therapy_";
Drugsbypatients.Current(idxtoreplace) = "_No Therapy_";
Drugsbypatients = removevars(Drugsbypatients,'patients');

% find (drug therapies, count) and display it (after doublons were discarded)
CFTR_modulator = categorical(unique(brDrugTherapy.DrugTherapyType));
druglisttable = table(CFTR_modulator);
for i = 1:size(druglisttable,1)
    druglisttable.Amount_prescribed(i) = sum(contains(Drugsbypatients.History,string(druglisttable.CFTR_modulator(i))));
end
druglisttable = sortrows(druglisttable,2,'descend');
fprintf('%i patients were prescribed at least one CFTR modulator:\n', length(patients));
disp(druglisttable)

% view patients with current drug therapy
currenttable =  groupcounts(Drugsbypatients,'Current');
currenttable = sortrows(currenttable,2,'descend');
fprintf('CFTR modulator status among the %i patients\n', size(Drugsbypatients,1));
disp(currenttable);

%%% find (patient CFTR modulator history, count, IDs) and display it
% group patients by their history
historytable =  groupcounts(Drugsbypatients,'History');

% get patient IDs for each drug mix
ID = string(zeros(size(historytable,1),1));
for i = 1:size(historytable,1)
    ID(i,1) = join(string(Drugsbypatients{ismember(Drugsbypatients.History,string(historytable.History(i))),'ID'}),', ');
end
historytable.IDs = ID;

% display results
fprintf('The %i patients have %i differents CFTR modulator history:\n', sum(historytable.GroupCount), size(historytable,1));
historytable = sortrows(historytable,2,'descend');
disp(historytable)

end