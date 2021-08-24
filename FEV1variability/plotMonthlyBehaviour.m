% plots monthly evolution of biomarkers

clc; clear; close all;
addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/smartcare/source_code/
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
study = 'BR';

[datamatfile, ~, ~] = getRawDataFilenamesForStudy(study);
[physdata, offset] = loadAndHarmoniseMeasVars(datamatfile, subfolder, study);


%% plot measurement values per month over years

% TODO % add number of people under tritherapy

year = 2019:1:2021;
month = 1:12;
%type = {'FEV1Recording'        }; feature_idx = 8;
%type = {'O2SaturationRecording'}; feature_idx = 10;
%type = {'RestingHRRecording'   }; feature_idx = 11;
%type = {'WellnessRecording'    }; feature_idx = 14;
type = [{'FEV1Recording'        }, {'O2SaturationRecording'        }, ...
    {'RestingHRRecording'   }, {'WellnessRecording'    }]; 
feature_idx = [8, 10, 11, 14];

f = figure('WindowState', 'fullscreen');
for i = 1:size(type,2)
    out = concatenateMonthlyData(physdata, year, month, feature_idx(i), type(i));
    labels = dateLabels(year, month);
    
    % don't plot over 1 months after latest date recorded
    [~, patientmasterdate, ~] = getLatestBreatheDatesForHosp('PAP');
    if str2num(patientmasterdate(1:4)) == year(end)
        last_month = size(year,2)*size(month,2) - 11 + str2num(patientmasterdate(5:6));
    end
    
    subplot(2,2,i);
    
    boxplot(out(:,1:last_month), 'Labels', labels(:,1:last_month))
    xlabel('Date')
    xtickangle(45)
    ylabel(type(i))
    title(['Monthly evolution of the ', type(i)])
    
    saveas(gcf,'../../PlotTristan/Monthly evolution of measures.png')
end


%% functions

function out = concatenateMonthlyData(physdata, year, month, feature_idx, type)
% concatenates the measurement values for a feature type, per month over
% the requested years
% input: data, year, months, feature, type, maxdate
% output: column vector with measurements values
    date = datevec(physdata.Date_TimeRecorded);
    
    % find indices for entries matching the conditions
    logical_array = compareLogicalArrays(date(:,1) == year, ...
        date(:,2) == month);
    logical_array = compareLogicalArrays(logical_array, ...
        ismember(physdata.RecordingType, type));

    % extracts FEV values for matched indices
    for i = 1:size(logical_array,2)
        filtered_val = physdata(find(logical_array(:,i)),feature_idx).(1);
        out(1:numel(filtered_val), i) = filtered_val;
    end
    % replace padded 0 by nan
    out(out == 0) = nan;
end

function logical_array = compareLogicalArrays(a, b)
% compares col by col multiple logical arrays
% input: nxm and nxp logical arrays
% output: (n,m*p) logical array
    logical_array = [];
    % for the m col of a
    for temp = a
        % perform an AND logic between the temp (nx1) and b (nxp)
        % thus adding p columns to the logical array
        logical_array = [logical_array, temp & b]; 
    end
end

function str = dateLabels(a,b);
% creates string labels based on year and months
    str = [];
    for i = 1:size(a,2)
        for j = 1:size(b,2)
            str = [str, sprintf("%i.%i ", b(j), a(i))];
        end
    end
end

% unused
function out = getMonthlyMeasures(data, year, month, feature)
    % input
    idx = 1:1:size(data,1);
    date = datevec(data.Date_TimeRecorded(idx));

    % find indices for entries matching the conditions
    matched_idx = find(ismember(date(idx,1), year) ... % year
        & ismember(date(idx,2), month) ... % month
        & ismember(data(idx,5).RecordingType, feature)); % feature

    % extracts FEV values for matched indices
    out = data(matched_idx,8).(1);
end