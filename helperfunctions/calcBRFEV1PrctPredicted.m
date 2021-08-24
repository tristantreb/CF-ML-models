% calculates FEV1 % predicted for each patient with values in Breathe data
% 
% formula: FEV1%Predicted = mean_of_all_FEV1_recordings / CalcPredictedFEV1
% 
% Input:
% ------
% ID, FEV1 measurement data
% ID, CalcPredictedFEV1 from cdPatient in the clinical data
% 
% Output: table with (ID, Value) with the percentage predicted
% -------

function FEV1prctpredicted = calcBRFEV1PrctPredicted(brphysdata,cdPatient)

fprintf('\nCalculating FEV1 percentage predicted');
% -> calculate PercentagePredicted based on brphysdata
FEV1breathe = getMeasureTable(brphysdata,'FEV1Recording','FEV');
% calculate true value based on mean
func = @(x) mean(x);
FEV1breathe = varfun(func,FEV1breathe,'GroupingVariables','ID','InputVariables','FEV1Recording');
% adds predicted FEV1
FEV1breathe = outerjoin(FEV1breathe, cdPatient, 'Type', 'Left', 'Keys', 'ID', 'RightVariables', {'CalcPredictedFEV1'});
% computes % predicted
FEV1breathe.PercentagePredicted = FEV1breathe.Fun_FEV1Recording ./ FEV1breathe.CalcPredictedFEV1 * 100;
FEV1breathe = renamevars(FEV1breathe, 'PercentagePredicted', 'Value');
FEV1prctpredicted = FEV1breathe(:,{'ID','Value'});
fprintf(' - Done!\n');

end

