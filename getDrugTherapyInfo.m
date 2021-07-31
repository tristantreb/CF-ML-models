function [Drugsbypatients] = getDrugTherapyInfo(brDrugTherapy) 
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

patients = unique(brDrugTherapy.ID);

%%% find (drug therapies, count) and display it
Drug_therapy = categorical(unique(brDrugTherapy.DrugTherapyType));
Count = countcats(categorical(brDrugTherapy.DrugTherapyType));
fprintf('Found %i drug therapies among %i patients:\n', sum(Count), length(patients))
disp(table(Drug_therapy,Count))

%%% find (drug mix, count, IDs) and display it

% group drug therapies by patients
list = string(zeros(length(patients),1));
for i = 1:length(patients)
    drugs_list = brDrugTherapy(brDrugTherapy.ID == patients(i),{'DrugTherapyStopDate', 'DrugTherapyType'});;
    list(i,1) = join(string(drugs_list.DrugTherapyType),', ');
    if ~isnat(drugs_list.DrugTherapyStopDate(end))
        list(i,1) = append(list(i,1),", _Therapy Stopped_");
    end
end

% (out var) table with patients and their associated drug mix
Drugsbypatients = table(patients,list);

% group patients by drug mix
Drugs_mix = categorical(unique(list));
Patient_count = countcats(categorical(list));

% get patient IDs for each drug mix
IDs = string(zeros(length(Drugs_mix),1));
for i = 1:length(Drugs_mix)
    IDs(i,1) = join(string(Drugsbypatients{ismember(Drugsbypatients.list,string(Drugs_mix(i))),'patients'}),', ');
end

% display results
fprintf('Found %i drug mix among %i patients:\n', length(Drugs_mix), length(patients))
disp(table(Drugs_mix,Patient_count,IDs))

end

