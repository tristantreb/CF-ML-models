% master script to run an recovery model

% Documentation Standard
%%%%%%%%%%%%%%%%%%%%%%%%
% "Purpose (1 line)"
%
% "Description"
% 
% Input:
% ------
% "var                  var definition (if not self-explicit)"
% "filename             what it contains (if not self-explicit)"
%
% Output:
% -------
% "variables, files, plots in the same way as input"

% CREATE DOCUMENTATION
%%%%%%%%%%%%%%%%%%%%%%
% Add m2html to path
%addpath /Users/tristan.trebaol/Documents/PDM/Project/m2html/
% run the following line, 'mfiles' should note be set to "./", else it

% raises a display issue
%m2html('mfiles','Code/','htmldir','Documentation','recursive','on', 'global','on');

% MASTER FILE TO RUN AN RECOVERY MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init;
% check source code root directory -> setBaseDir

% CLINICAL, ID and StudyNumber (gmail)
% now automaticaly loads most recent file
% check folder MatlabSavedVariables exists
%loadbreatheclinicalREDCapdata;

% MAGIC BULLET, guid (Globally Unique Identifier)
% loads raw data 
%breatheAzureAccess;

%check date input -> getLatestBreatheMeasDate
%loadbreathemeasurementdata; %datademographics

% creates treatment lists - merge concurrent treatments
%checkBreatheTreatmentsVsSelectMeasuresPost;
%gap in treatments

% confront clinical vs home measures
% weight
%checkBreatheClinicalVsHomeWeightMeasures;
% fev1
%checkBreatheClinicalVsHomeFEV1Measures;

%createAlignmentModelInputsRecovery;
% study number
% gap in treatments

%findDataAnomaliesRecovery;
% study number
% gap in treatments

% missing files
% BR_LabelledInterventions_gap10consensus.mat - given
% BRelectivetreatments_gap10.xlsx - given
%runAlignmentModelEMMCScriptRecovery;