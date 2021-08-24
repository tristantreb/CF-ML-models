% master script to run an alignment model for APEs (prior to redcap)

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

% run the following line, 'mfiles' should note be set to "./", else it
% raises a display issue
% addpath /Users/tristan.trebaol/Documents/PDM/Project/m2html/
% m2html('mfiles','Code/','htmldir','Documentation/','recursive','on', 'global','on','todo','on');


% MASTER FILE TO RUN AN ALIGNMENT MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/
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
%checkBreatheTreatmentsVsSelectMeasuresPrior;
%checkIVTreatmentsVsMeasuresPriorNew; used before 7.05.2021
%gap in treatments

createAlignmentModelInputs;
% study number
% gap in treatments

findDataAnomalies;
% study number
% gap in treatments

% missing files
% BR_LabelledInterventions_gap10consensus.mat - given
% BRelectivetreatments_gap10.xlsx - given
runAlignmentModelEMMCScript;