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

addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/
% check source code root directory -> setBaseDir

% CLINICAL, ID and StudyNumber (gmail)
% check dates input -> getLatestBreatheDatesForHosp
% check folder MatlabSavedVariables exists
%loadbreatheclinicaldatanew;

% MAGIC BULLET, guid (Globally Unique Identifier)
% loads raw data 
%breatheAzureAccess;

%check date input -> getLatestBreatheMeasDate
%loadbreathedata; %datademographics

% creates treatment lists - merge concurrent treatments
%checkBreatheTreatmentsVsSelectMeasuresPost;
%gap in treatments

%createAlignmentModelInputsRecovery;
% study number
% gap in treatments

%findDataAnomaliesRecovery;
% study number
% gap in treatments

% missing files
% BR_LabelledInterventions_gap10consensus.mat - given
% BRelectivetreatments_gap10.xlsx - given
runAlignmentModelEMMCScriptRecovery;