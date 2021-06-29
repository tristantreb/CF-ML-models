clear;
addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/smartcare/;
addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/pdm/;
addpath /Users/tristan.trebaol/Documents/PDM/Project/Code/pdm/recovery/;
basedir = setBaseDir();
subfolder = 'MatlabSavedVariables';
study = 'BR';
plotfolder = '/Users/tristan.trebaol/Documents/PDM/Project/PlotsRecovery/';

chosentreatgap = 10;
data_window = 40;

% colors
cvcol   = [0.94  0.52  0.15];
admcol  = [0.694 0.627 0.78]; 
ivcol   = [1     0     0   ];
oralcol = [1     0.85  0   ];
trplcol = [0     1     0   ];
drugcol = [0     0.8   0.6 ];

% manually add triple therapy start for 501
% manually add predicted FEV1 to 3.69 for 616