function [out] = getMeasureTable(brphysdata,recording_type,category)
% returns a table with all data about the specified recording type
%
% Input: raw data, recording_type, category of the recording type
% ------
% there are 17 options for the recording type:
%     'CalorieRecording'
%     'CoughRecording'
%     'FEF2575Recording'
%     'FEV075Recording'
%     'FEV1DivFEV6Recording'
%     'FEV1Recording'
%     'FEV6Recording'
%     'HasColdOrFluRecording'
%     'HasHayFeverRecording'
%     'MinsAsleepRecording'
%     'MinsAwakeRecording'
%     'O2SaturationRecording'
%     'PulseRateRecording'
%     'RestingHRRecording'
%     'TemperatureRecording'
%     'WeightRecording'
%     'WellnessRecording'
%
% Output: table with all information
% -------

colnames = {'SmartCareID','ScaledDateNum','DateNum',category,'UserName','CaptureType','Date_TimeRecorded'};
out = brphysdata(ismember(brphysdata.RecordingType, recording_type),colnames);
end

